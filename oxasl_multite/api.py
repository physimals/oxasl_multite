"""
OXASL plugin for processing multiphase ASL data

Copyright (c) 2019 Univerisity of Oxford
"""
import math
import numpy as np

from fsl.wrappers import LOAD
from fsl.data.image import Image

from oxasl import basil
from oxasl.options import OptionCategory, IgnorableOptionGroup
from oxasl.reporting import Report
from oxasl.wrappers import fabber, mvntool

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

def _base_fabber_options(wsp, asldata):
    options = {
        "method" : "vb",
        "noise" : "white",
        "model" : "asl_multite",
        "data" : asldata,
        "mask" : wsp.rois.mask,
        "ti" : list(asldata.tis),
        "tau" : list(asldata.taus),
        "repeats" : asldata.rpts[0], # We have already checked repeats are fixed
        "save-mean" : True,      
        "max-iterations": 30,
    }

    if wsp.bat is None:
        wsp.bat = 1.3 if wsp.casl else 0.7
    if wsp.batsd is None:
        wsp.batsd = 1.0 if len(asldata.tis) > 1 else 0.1

    for opt in ("bat", "batsd", "t1", "t1b"):
        val = wsp.ifnone(opt, None)
        if val is not None:
            options[opt] = val

    return options

def _multite_fabber_options(wsp, asldata):
    options = _base_fabber_options(wsp, asldata)
    options.update({
        "model" : "asl_multite",
        "te" : list(wsp.asldata.tes),
        "infertexch" : True,
        "save-model-fit" : True,        
        "max-iterations": 30,
    })
    # Additional user-specified multiphase fitting options override the above
    options.update(wsp.ifnone("multite_options", {}))
    return options

def _aslrest_fabber_options(wsp, asldata):
    options = _base_fabber_options(wsp, asldata)
    options.update({
        "model" : "aslrest",
        "casl" : True,
        "inctiss" : True,
        "incbat" : True,
        "infertiss" : True,
        "inferbat" : True,
        "save-std" : True,
    })
    return options

def fit_init(wsp):
    """
    Do an initial fit on ftiss and delttiss using the aslrest model
    """

    # Estimate the FTISS/DELTISS using the first TE value. 
    # This is easiest to extract if the TEs are slowest varying
    # Note that it would be better to use all of the TEs but we need 
    # to correct the output for a single TE so that wouldn't work.
    # What *would* work is applying the T2 correction to the input
    # data - will consider that for the future. 
    data_multite = wsp.asldata.diff().reorder(out_order="tre").data
    ntes = wsp.asldata.ntes
    nvols_mean = int(wsp.asldata.nvols/ntes)
    data_mean = np.zeros(list(data_multite.shape[:3]) + [nvols_mean])
    for idx in range(1): # FIXME first TE only?
        data_mean += data_multite[..., idx*nvols_mean:(idx+1)*nvols_mean]

    wsp.asldata_mean = wsp.asldata.derived(image=data_mean, name="asldata", 
                                           iaf="diff", order="tr", tes=[0])

    # Run ASLREST on the mean data to generate initial estimates for CBF and ATT
    # Note that the multi-TE model has an additional T2 correction factor of exp(-te/T2) which ASLREST lacks.
    # Since we are using the shortest TE we apply this correction to the FTISS estimates
    options = _aslrest_fabber_options(wsp, wsp.asldata_mean)
    result = _run_fabber(wsp.sub("aslrest"), options, "Running Fabber using standard ASL model for CBF/ATT initialization")
    t2_corr_factor = math.exp(-wsp.asldata.tes[0] / wsp.t2)
    wsp.log.write("  - Using T2 correction factor on ASLREST output: %f\n" % t2_corr_factor)
    wsp.aslrest.mean_ftiss_t2corr = Image(wsp.aslrest.mean_ftiss.data / t2_corr_factor, header=wsp.aslrest.mean_ftiss.header)
    wsp.aslrest.var_ftiss = Image(np.square(wsp.aslrest.std_ftiss.data / t2_corr_factor), header=wsp.aslrest.std_ftiss.header)
    wsp.aslrest.var_delttiss = Image(np.square(wsp.aslrest.std_delttiss.data), header=wsp.aslrest.std_delttiss.header)

    # Run the multi-TE model for 1 iteration to get an MVN in the correct format
    options = _multite_fabber_options(wsp, wsp.asldata)
    options.update({"save-mvn" : True, "max-iterations" : 1})
    result = _run_fabber(wsp.sub("mvncreate"), options, "Running Fabber for 1 iteration on multi-TE model to generate initialization MVN")

    # Merge the CBF and ATT estimates from the ASLREST run into the output MVN to generate an initialization MVN
    # for the final multi-TE fit.
    wsp.log.write("  - Merging CBF and ATT estimates into the MVN to initialize multi-TE fit\n")
    wsp.init_mvn = mvntool(wsp.mvncreate.finalMVN, 1, output=LOAD, mask=wsp.rois.mask, write=True, valim=wsp.aslrest.mean_ftiss_t2corr, varim=wsp.aslrest.var_ftiss, log=wsp.fsllog)["output"]
    wsp.init_mvn = mvntool(wsp.init_mvn, 2, output=LOAD, mask=wsp.rois.mask, write=True, valim=wsp.aslrest.mean_delttiss, varim=wsp.aslrest.var_delttiss, log=wsp.fsllog)["output"]

def fit_multite(wsp):
    """
    """
    wsp.log.write("\nPerforming multi-TE model fitting:\n")
    if wsp.asldata.is_var_repeats():
        raise ValueError("Multi-TE ASL data with variable repeats not currently supported")

    # Make sure repeats are the slowest varying as this is what the model expects. Similarly
    # make sure varying TEs are always within each TI
    wsp.asldata = wsp.asldata.diff().reorder(out_order="etr")
    options = _multite_fabber_options(wsp, wsp.asldata)

    if wsp.multite_init:
        wsp.sub("init")
        fit_init(wsp.init)
        options["continue-from-mvn"] = wsp.init.init_mvn

    result = _run_fabber(wsp.multite.sub("finalstep"), options, "Running Fabber using multi-TE model")
    wsp.log.write("\nDONE multi-TE model fitting\n")

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
