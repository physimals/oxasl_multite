"""
OXASL plugin for processing multiphase ASL data

Copyright (c) 2019 Univerisity of Oxford
"""
import numpy as np

from fsl.wrappers import LOAD
from fsl.data.image import Image

from oxasl import basil
from oxasl.options import OptionCategory, IgnorableOptionGroup
from oxasl.reporting import Report
from oxasl.wrappers import fabber

from ._version import __version__

def _run_fabber(wsp, options, desc):
    """
    Run Fabber and write the output to a workspace
    """
    wsp.log.write("  - %s     " % desc)
    result = fabber(options, output=LOAD, progress_log=wsp.log, log=wsp.fsllog)
    wsp.log.write(" - DONE\n")

    for key, value in result.items():
        setattr(wsp, key, value)

    if result["logfile"] is not None and wsp.savedir is not None:
        wsp.set_item("logfile", result["logfile"], save_fn=str)
    return result

def _fabber_options(wsp):
    """
    :return: General Fabber options for multi-TE decoding
    """

    # General options. note that the phase is always PSP number 1
    options = {
        "method" : "vb",
        "noise" : "white",
        "model" : "asl_multite",
        "data" : wsp.asldata,
        "mask" : wsp.rois.mask,
        "ti" : list(wsp.asldata.tis),
        "te" : list(wsp.asldata.tes),
        "tau" : list(wsp.asldata.taus),
        "repeats" : wsp.asldata.rpts[0], # We have already checked repeats are fixed
        "infertexch" : True,
        "save-mean" : True,
        "save-model-fit" : True,        
        "max-iterations": 30,
    }

    for opt in ("t1", "t1b", "t2", "t2b"):
        val = wsp.ifnone(opt, None)
        if val is not None:
            options[opt] = val

    # Additional user-specified multiphase fitting options override the above
    options.update(wsp.ifnone("multite_options", {}))

    return options

def fit_multite(wsp):
    """
    """
    wsp.log.write("\nPerforming multi-TE model fitting:\n")
    if wsp.asldata.is_var_repeats():
        raise ValueError("Multi-TE ASL data with variable repeats not currently supported")

    # Make sure repeats are the slowest varying as this is what the model expects. Similarly
    # make sure varying TEs are always within each TI
    wsp.asldata = wsp.asldata.diff().reorder(out_order="etr")

    if wsp.multite_init:
        wsp.sub("init")
        fit_init(wsp.init)
        init_mvn = wsp.init.initial_mvn
    else:
        init_mvn = None

    # Get the Fabber options
    options = _fabber_options(wsp.multite)
    print(options)
    result = _run_fabber(wsp.multite.sub("finalstep"), options, "Running Fabber using multi-TE model")

    wsp.log.write("\nDONE multi-TE decoding\n")

def model_multite(wsp):
    """
    Do modelling on multi-TE ASL data

    :param wsp: Workspace object

    Required workspace attributes
    -----------------------------

      - ``asldata`` - ASLImage containing multi-TE data

    Optional workspace attributes
    -----------------------------

    See ``MultiTEOptions`` for other options

    Workspace attributes updated
    ----------------------------

      - ``multite``    - Sub-workspace containing multi-TE decoding output
      - ``output``     - Sub workspace containing native/structural/standard space
                         parameter maps
    """
    wsp.sub("multite")
    fit_multite(wsp.multite)

    # Write output
    wsp.sub("output")

    from oxasl import oxford_asl
    oxford_asl.output_native(wsp.output, wsp.multite)

    # Re-do registration using PWI map.
    oxford_asl.redo_reg(wsp, wsp.output.native.perfusion)

    # Write output in transformed spaces
    oxford_asl.output_trans(wsp.output)

    wsp.log.write("\nDONE processing\n")

class MultiTEOptions(OptionCategory):
    """
    OptionCategory which contains options for preprocessing multi-TE ASL data
    """
    def __init__(self, **kwargs):
        OptionCategory.__init__(self, "oxasl_multite", **kwargs)

    def groups(self, parser):
        groups = []
        group = IgnorableOptionGroup(parser, "Multi-TE Options", ignore=self.ignore)
        group.add_option("--multite-init", help="Initialize perfusion and transit time using fit on restring state ASL model", action="store_true", default=False)
        group.add_option("--multite-options", help="File containing additional options for multiphase fitting step", type="optfile")
        groups.append(group)
        return groups
