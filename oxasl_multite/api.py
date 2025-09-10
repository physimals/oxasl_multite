"""
OXASL plugin for processing multiphase ASL data

Copyright (c) 2019 Univerisity of Oxford
"""
import math
import numpy as np

from fsl.wrappers import LOAD
from fsl.data.image import Image

from oxasl import basil
from oxasl.options import OptionCategory, OptionGroup
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
        "nte" : len(wsp.asldata.tes),
        "infertexch" : True,
        "save-std" : True,
        "save-model-fit" : True,
        "save-residuals" : wsp.ifnone("output_residuals", False),
        "max-iterations": 30,
        "t2" : float(wsp.ifnone("t2", 50)) / 1000
    })

    if wsp.spatial:
        options.update({
            "PSP_byname1" : "ftiss",
            "PSP_byname1_type" : "M",
            "method" : "spatialvb",
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

def init_t2(wsp):
    """
    Initialize the T2 value by fitting the T2 decay part of the signal

    We do not use Fabber for this (although it would be possible). Instead
    we do a simple voxel-by-voxel least squares fit to a T2 decay model using
    a subset of voxels with the strongest signals. We take the median of the T2 
    value as our estimate since the mean can be affected by extreme values resulting
    from fitting problems.
    """
    def t2model(tes, t2, *s0):
        ntes = len(tes)/len(s0)
        s0 = np.repeat(np.array(s0), ntes)
        return s0 * np.exp(-1000*tes/t2)

    wsp.log.write("  - Initializing T2 value by fit to T2 decay model\n")
    
    # Estimate T2 and a T2 corrected signal
    wsp.data_multite = wsp.asldata.diff().reorder(out_order="etr")
    data_multite = wsp.data_multite.data
    volshape = list(data_multite.shape[:3])
    
    # Identify unmasked voxels with strongest signals
    diffdata = data_multite.max(axis=-1) - data_multite.min(axis=-1)
    thresh = np.percentile(diffdata, wsp.ifnone("multite_t2_init_percentile", 95))
    wsp.log.write("  - Including unmasked voxels with signal > %f\n" % thresh)
    mask_thresh = diffdata > thresh
    roidata = np.logical_and(wsp.rois.mask.data > 0, diffdata > thresh)

    data_multite_roi = data_multite[roidata]
    nvoxels_roi = data_multite_roi.shape[0]
    nvols = wsp.asldata.nvols
    ntes = wsp.asldata.ntes
    nsigs = int(nvols / ntes)
    tes = np.array(wsp.asldata.tes * nsigs)
    wsp.log.write("  - %i TEs, %i volumes, %i signals, %i voxels\n" % (ntes, nvols, nsigs, nvoxels_roi))

    # Go through each voxel and fit the T2 decay model for S0 and T2
    t2_roi = np.zeros((nvoxels_roi, ), dtype=np.float32)
    sig_roi = np.zeros((nvoxels_roi, nsigs), dtype=np.float32)
    for voxel_idx, data_multite_voxel in enumerate(data_multite_roi):
        try:
            # Initialize signal from maximum of the data at each time point
            sig_init = [max(data_multite_voxel[sig_idx*ntes:(sig_idx+1)*ntes]) for sig_idx in range(nsigs)]
            param_init=[wsp.t2, ] + sig_init
            from scipy.optimize import curve_fit
            popt, pcov = curve_fit(t2model, tes, data_multite_voxel, p0=param_init)
            t2_roi[voxel_idx] = popt[0]
            sig_roi[voxel_idx, :] = popt[1:]
        except Exception as exc:
            wsp.log.write("  - WARNING: fit failed for voxel: %i\n" % voxel_idx)
    wsp.t2 = np.median(t2_roi)
    wsp.log.write("  - Median T2: %f ms\n" % wsp.t2)

def fit_init(wsp):
    """
    Do an initial fit on ftiss and delttiss using the aslrest model

    The first stage of this is to apply a T2 correction to the multi-TE data
    and take the mean across TEs, since the ASLREST model contains no T2
    correction.

    The resulting ASL data is put through the basic ASLREST model. We then run a single
    iteration of the multi-TE model to generate an MVN, and insert the ftiss and delttiss
    output from ASLREST into the MVN. This is then used to initialize the multi-TE run.
    """
    wsp.log.write("  - Preparing initialization of perfusion from resting state model\n")
    wsp.data_multite = wsp.asldata.diff().reorder(out_order="etr")
    data_multite = wsp.data_multite.data
    tes = wsp.asldata.tes
    nvols_mean = int(wsp.asldata.nvols/len(tes))

    # Do the T2 correction and take the mean across TEs
    data_mean = np.zeros(list(data_multite.shape[:3]) + [nvols_mean])
    for idx, te in enumerate(tes):
        t2_corr_factor = math.exp(1000 * te / wsp.t2)
        wsp.log.write("  - Using T2 correction factor for TE=%f: %f\n" % (te, t2_corr_factor))
        data_mean += data_multite[..., idx::wsp.asldata.ntes] * t2_corr_factor
    wsp.asldata_mean = wsp.asldata.derived(image=data_mean/len(tes), name="asldata", 
                                           iaf="diff", order="tr", tes=[0])

    # Run ASLREST on the mean data to generate initial estimates for CBF and ATT
    options = _aslrest_fabber_options(wsp, wsp.asldata_mean)
    result = _run_fabber(wsp.sub("aslrest"), options, "Running Fabber using standard ASL model for CBF/ATT initialization")
    wsp.aslrest.var_ftiss = Image(np.square(wsp.aslrest.std_ftiss.data), header=wsp.aslrest.std_ftiss.header)
    wsp.aslrest.var_delttiss = Image(np.square(wsp.aslrest.std_delttiss.data), header=wsp.aslrest.std_delttiss.header)

    # Run the multi-TE model for 1 iteration to get an MVN in the correct format
    options = _multite_fabber_options(wsp, wsp.asldata)
    options.update({"save-mvn" : True, "max-iterations" : 1})
    result = _run_fabber(wsp.sub("mvncreate"), options, "Running Fabber for 1 iteration on multi-TE model to generate initialization MVN")

    # Merge the CBF and ATT estimates from the ASLREST run into the output MVN to generate an initialization MVN
    # for the final multi-TE fit.
    wsp.log.write("  - Merging CBF and ATT estimates into the MVN to initialize multi-TE fit\n")
    wsp.init_mvn = mvntool(wsp.mvncreate.finalMVN, 1, output=LOAD, mask=wsp.rois.mask, write=True, valim=wsp.aslrest.mean_ftiss, varim=wsp.aslrest.var_ftiss, log=wsp.fsllog)["output"]
    wsp.init_mvn = mvntool(wsp.init_mvn, 2, output=LOAD, mask=wsp.rois.mask, write=True, valim=wsp.aslrest.mean_delttiss, varim=wsp.aslrest.var_delttiss, log=wsp.fsllog)["output"]

def fit_multite(wsp):
    """
    Run model fitting on multi-TE data
    """
    wsp.log.write("\nPerforming multi-TE model fitting:\n")
    if wsp.asldata.is_var_repeats():
        raise ValueError("Multi-TE ASL data with variable repeats not currently supported")

    # Make sure repeats are the slowest varying as this is what the model expects. Similarly
    # make sure varying TEs are always within each TI
    wsp.asldata = wsp.asldata.diff().reorder(out_order="etr")
    options = _multite_fabber_options(wsp, wsp.asldata)

    if wsp.multite_init_t2:
        wsp.sub("init_t2")
        init_t2(wsp.init_t2)
        wsp.t2 = wsp.init_t2.t2

    if wsp.multite_init:
        wsp.sub("init")
        fit_init(wsp.init)
        options["continue-from-mvn"] = wsp.init.init_mvn

    result = _run_fabber(wsp.multite.sub("finalstep"), options, "Running Fabber using multi-TE model")
    wsp.log.write("\nDONE multi-TE model fitting\n")

def run(wsp):
    fit_multite(wsp.sub("multite"))
    wsp.quantify_wsps.append("multite")

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
        group = OptionGroup(parser, "Multi-TE Options")
        group.add_option("--multite-init-t2", help="Initialize T2 value", action="store_true", default=False)
        group.add_option("--multite-init", help="Initialize perfusion and transit time using fit on restring state ASL model", action="store_true", default=False)
        group.add_option("--multite-options", help="File containing additional options for multiphase fitting step", type="optfile")
        groups.append(group)
        return groups
