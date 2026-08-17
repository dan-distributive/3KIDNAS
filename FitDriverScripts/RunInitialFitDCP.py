"""
RunInitialFitDCP.py
====================

Runs the ONE initial/anchor fit for a galaxy (Fortran's equivalent:
Programs/SingleGalaxyFitter run with BSSwitch=0, via RunWRKP.py) entirely
through the JS/DCP pipeline instead --
js/app/bootstrap-realization-launcher.js's runInitialFit(). Mirrors
RunBootstrapsDCP.py's shape (payload builder / dispatch / result adapter),
but is a SEPARATE module, not a shared function -- deliberately does not
modify RunBootstrapsDCP.py (itself still marked "NOT YET RUN END-TO-END,
review carefully" in its own header) to avoid any regression risk to the
already-working bootstrap-realization dispatch.

Called from FullSingleGalaxyFit.GalaxyFit() when GalaxyDict['UseDCP']==1,
INSTEAD of RW.RunWRKP(GeneralDict, GalaxyDict, 0) + RWF.ReadWRKPOutputFile
-- sets GalaxyDict['BestFitModel'] directly, no Fortran, no text-file
round-trip. UseDCP==0 (fortran-local) is completely untouched by this file.

Unlike a bootstrap realization, the initial fit has no pre-existing fixture
to read cube/beam geometry from (that fixture is normally PRODUCED by the
initial fit) -- runInitialFit() derives pixel/channel/beam geometry itself,
directly from the raw FITS bytes, so this payload builder only needs to
supply the raw cube+mask bytes, the static PA/Inc estimate, and fitting-
options config values -- see runInitialFit()'s own header comment in
bootstrap-realization-launcher.js for the full precision rationale.
"""

import os
import json
import subprocess
import base64
import warnings
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.wcs import FITSFixedWarning


# Fortran PipelineGlobals.f hardcoded parameter constants -- not read from
# any config file (same constants RunBootstrapsDCP.py's V_ROT_LIMS/SIZE_LIMS
# replicate; LIN_SD_LIMS_MSOL is NOT needed here -- see below, that one's
# cube-geometry-dependent derivation now happens in JS instead).
V_ROT_LIMS = [0.0, 400.0]
SIZE_LIMS = [2.0, 10.0]


# ---------------------------------------------------------------------------
# 1. Fitting-options file parsing -- full sequential read order, replaying
#    src/Inputs/FittingOptionsInputs.f's FittingOptionsIn() over the SAME
#    options-file lines Python already loads (GeneralDict
#    ['MainWRKPOptionsLines'], via RunWRKP.LoadDefaultWRKPFiles). Continues
#    further than RunBootstrapsDCP.ParseFittingOptionsExtras does (that one
#    stops after the ring-grid block, since a bootstrap realization gets
#    sigmaLengths/nRingsPerBeam/constParams/fixedParams from the fixture
#    instead) -- the initial fit has no fixture, so needs all of them read
#    directly from the options file, matching FittingOptionsInputs.f:1-220
#    exactly, confirmed by direct read of that file.
# ---------------------------------------------------------------------------
def ParseFittingOptionsForInitialFit(optionsLines):
    it = iter(optionsLines)

    def read_value():
        next(it)          # label/comment line
        return next(it).strip()

    read_value()  # core_code_switch -- walked, not used (JS always uses GalaxyFit_Simple)
    likelihood_switch = int(read_value())
    read_value()  # param_interface_switch -- walked, not used (JS always uses GeneralizedParamVectorToTiltedRing)
    moment_map_switch = int(read_value())

    if moment_map_switch == 3:
        next(it); next(it); next(it)             # label, SecondaryCubeFile, SecondaryMaskFile
    elif moment_map_switch == 2:
        next(it); next(it); next(it); next(it); next(it)  # label, 3 MapNames, ProfileName

    center_source = int(read_value())
    read_value()  # shape_source -- walked, not used (initial fit's geometry is always the static catalogue estimate)
    read_value()  # size_source -- walked, not used
    read_value()  # vsys_source -- walked, not used
    sd_switch = int(read_value())
    cmode = int(read_value())
    cloud_base_surf_dens = float(read_value())
    sigma_lengths = float(read_value())
    noise_sigma_lim = float(read_value())
    n_targ_rings_raw = int(read_value())

    rad_grid_arcsec = None
    if n_targ_rings_raw > 0:
        next(it)  # label
        rad_str_line = next(it)
        rad_grid_arcsec = [float(s) for s in rad_str_line.split()]

    n_rings_per_beam = int(read_value())

    const_params = [False] * 13
    fixed_params = [False] * 13
    for i in range(13):
        next(it)  # label line
        toks = next(it).split()
        const_params[i] = toks[0].strip().upper().startswith('T')
        fixed_params[i] = toks[1].strip().upper().startswith('T')

    return {
        'centerSource': center_source,
        'sdSwitch': sd_switch,
        'cmode': cmode,
        'cloudBaseSurfDens': cloud_base_surf_dens,
        'sigmaLengths': sigma_lengths,
        'noiseSigmaLim': noise_sigma_lim,
        'nTargRings': n_targ_rings_raw,
        'radGridArcsec': rad_grid_arcsec,
        'nRingsPerBeam': n_rings_per_beam,
        'constParams': const_params,
        'fixedParams': fixed_params,
        'likelihoodSwitch': likelihood_switch,
    }


# ---------------------------------------------------------------------------
# 2. Build the payload for the single initial-fit job.
# ---------------------------------------------------------------------------
def BuildInitialFitPayload(GeneralDict, GalaxyDict):
    extras = ParseFittingOptionsForInitialFit(GeneralDict['MainWRKPOptionsLines'])

    with open(GalaxyDict['CubeNameU'], 'rb') as f:
        observed_raw_fits_b64 = base64.b64encode(f.read()).decode('ascii')
    with open(GalaxyDict['MaskNameU'], 'rb') as f:
        mask_raw_fits_b64 = base64.b64encode(f.read()).decode('ascii')

    return {
        'observedCubeRawFitsB64': observed_raw_fits_b64,
        'observedMaskRawFitsB64': mask_raw_fits_b64,
        # Static geometry -- matches Fortran's own initial fit exactly:
        # MakeCatalogue()/PFlags%CatFlag just reads these two numbers off a
        # 2-line text file (SoFiA_Driver.WriteSoFiACatFileForWRKP) -- no
        # live SoFiA run for the initial fit, in Fortran either.
        'paEstDeg': float(GalaxyDict['PA_EstimateU']),
        'incEstDeg': float(GalaxyDict['Inc_EstimateU']),
        'centerSource': extras['centerSource'],
        'nRingsPerBeam': extras['nRingsPerBeam'],
        'nTargRings': extras['nTargRings'],
        'radGridArcsec': extras['radGridArcsec'],
        'sdSwitch': extras['sdSwitch'],
        'constParams': extras['constParams'],
        'fixedParams': extras['fixedParams'],
        'vRotLims': V_ROT_LIMS,
        'sizeLims': SIZE_LIMS,
        'noiseSigmaLim': extras['noiseSigmaLim'],
        'sigmaLengths': extras['sigmaLengths'],
        'cmode': extras['cmode'],
        'cloudBaseSurfDens': extras['cloudBaseSurfDens'],
        'fitIdum': ReadMainInputSeed(GeneralDict['MainWRKPInputLines']),
        # galaxyFit_Simple hardcodes ftol=0.005 for pass 1 regardless of what
        # is passed in (see GalaxyFit.js / src/GalaxyAnalysis/GalaxyFit.f:93)
        # -- this value only ever matters as the pass-2 divisor's input
        # (ftol/5.), so 0.005 here reproduces the real starting point exactly.
        'ftol': 0.005,
        'likelihoodSwitch': extras['likelihoodSwitch'],
    }


# Random-seed line in the main WRKP input template (SingleFitInput_Base.in
# line 12, a FIXED pipeline default -- e.g. "-3" -- NOT derived from
# GalaxyDict['BootstrapSeed'], which is a separate, resample-only seed read
# by a different part of the config (BootstrapRuntimeInputs.f). Read
# directly off the same GeneralDict['MainWRKPInputLines'] Fortran itself
# reads idum from (RunWRKP.WriteWRKPMainFile never overwrites this line),
# rather than assumed/hardcoded, so a differently-configured
# WRKP_GeneralMainIn template is still read correctly.
def ReadMainInputSeed(mainInputLines):
    return int(mainInputLines[11].strip())


# ---------------------------------------------------------------------------
# 3. Launch the single DCP job.
# ---------------------------------------------------------------------------
def LaunchInitialFitJob(GeneralDict, payload, objName=''):
    workdir = GeneralDict['DCPJobDir']
    launcher = GeneralDict['DCPRealizationLauncher']
    os.makedirs(workdir, exist_ok=True)

    payload_path = os.path.join(workdir, 'initial_fit_payload.json')
    results_path = os.path.join(workdir, 'dcp_initial_fit_result.json')

    if os.path.exists(results_path):
        os.remove(results_path)

    with open(payload_path, 'w') as f:
        json.dump(payload, f)

    cmd = [
        'node', launcher,
        '--initialFit',
        '--payload', payload_path,
        '--out', results_path,
        '--jobName', '3KIDNAS-' + objName,
    ]
    # Same local/real-DCP dual-mode selection as RunBootstrapsDCP.
    # LaunchDCPRealizationJob -- see that function's own comment for why
    # DCP_FORCE_LOCAL exists.
    apiKey = os.environ.get('DCP_API_KEY')
    forceLocal = os.environ.get('DCP_FORCE_LOCAL') == '1'
    if apiKey and not forceLocal:
        cmd += ['--apiKey', apiKey]
    else:
        if not apiKey:
            print("DCP_API_KEY not set -- falling back to --local (direct in-process execution, no dispatch)")
        # --local expects a following value (bootstrap-realization-launcher.
        # js's shared `arg()` parser takes argv[i+1] unconditionally) -- the
        # initial-fit CLI branch never reads that value, just --local's
        # presence, but a bare trailing --local with nothing after it parses
        # as undefined, which arg() treats as "flag absent".
        cmd += ['--local', '1']
    computeGroups = os.environ.get('DCP_COMPUTE_GROUPS')
    if computeGroups:
        cmd += ['--computeGroups', computeGroups]
    slicePrice = os.environ.get('DCP_SLICE_PRICE')
    if slicePrice:
        cmd += ['--slicePrice', slicePrice]

    r = subprocess.run(cmd, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError("DCP initial-fit job failed (returncode %d)" % r.returncode)

    with open(results_path, 'r') as f:
        return json.load(f)


def GetCubeWCS(cube_path):
    hdr = fits.getheader(cube_path)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        w = wcs.WCS(hdr)
    return w


# ---------------------------------------------------------------------------
# 4. Adapter: worker result -> GalaxyDict['BestFitModel']-shaped dict,
#    matching ReadWRKPFit.LoadBestFitModelFile's contract (same keys
#    downstream code -- GeometryCorrection/ExtendedSDProfileCalc/
#    ExtractScalingParams/Bootstrap_Outputs -- already reads regardless of
#    which leg produced the fit). Also persists the two files that contract
#    requires to exist on disk: the model-cube FITS (SavePVDiagrams/
#    MakeBootstrapModelPlot fits.open() it directly) and diskfit_fixture.
#    json (RunBootstrapsDCP.BuildRealizationPayload reads it for the
#    SUBSEQUENT bootstrap loop).
# ---------------------------------------------------------------------------
def AdaptInitialFitResult(result, cube_wcs, GeneralDict, GalaxyDict):
    if result.get('sofiaFailed') or not result.get('FITAchieved', False):
        return {'FITAchieved': False}

    xc = np.array(result['XCENTER'], dtype=float)
    yc = np.array(result['YCENTER'], dtype=float)
    try:
        ra, dec = cube_wcs.celestial.all_pix2world(xc, yc, 0)
    except Exception:
        ra = np.full_like(xc, np.nan)
        dec = np.full_like(yc, np.nan)

    n = len(xc)
    zeros = np.zeros(n)
    model = {
        'FITAchieved': True,
        'CHI2': result['CHI2'],
        'XCENTER': list(xc), 'XCENTER_ERR': zeros,
        'YCENTER': list(yc), 'YCENTER_ERR': zeros,
        'INCLINATION': np.array(result['INCLINATION'], dtype=float), 'INCLINATION_ERR': zeros,
        'POSITIONANGLE': np.array(result['POSITIONANGLE'], dtype=float), 'POSITIONANGLE_ERR': zeros,
        'VSYS': np.array(result['VSYS'], dtype=float), 'VSYS_ERR': zeros,
        'VROT': np.array(result['VROT'], dtype=float), 'VROT_ERR': zeros,
        'VDISP': np.array(result['VDISP'], dtype=float), 'VDISP_ERR': zeros,
        'SURFDENS': np.array(result['SURFDENS'], dtype=float), 'SURFDENS_ERR': zeros,
        'SURFDENS_FACEON': np.array(result['SURFDENS_FACEON'], dtype=float), 'SURFDENS_FACEON_ERR': zeros,
        'RA': list(ra), 'RA_ERR': zeros,
        'DEC': list(dec), 'DEC_ERR': zeros,
        'R': np.array(result['R'], dtype=float),
        'R_SD': np.array(result['R_SD'], dtype=float),
        'RMS': np.full(1, result['RMS']),
        'SN_Integrated': np.full(1, result['SN_Integrated']),
        'SN_Peak': np.full(1, result['SN_Peak']),
        'SN_Avg': np.full(1, result['SN_Avg']),
        'SN_Median': np.full(1, result['SN_Median']),
    }
    # All the *_ERR fields above are placeholders (0) -- confirmed dead:
    # Bootstrap_Error_Analysis.EstimateUncertaintiesFromBootstraps
    # unconditionally overwrites every one of them from the bootstrap
    # ensemble before anything reads them (BootsrapErrors_ForParam).

    TargFolder = GalaxyDict['TargFolderU'] + "/" + GalaxyDict['ObjNameU'] + "/"
    os.makedirs(TargFolder, exist_ok=True)
    model_cube_path = TargFolder + GalaxyDict['ObjNameU'] + "_AverageModel_v1.fits"
    with open(model_cube_path, 'wb') as f:
        f.write(base64.b64decode(result['modelCubeFitsB64']))
    model['ModelCube'] = model_cube_path

    # diskfit_fixture.json -- same path GeneralDict['FixtureFile'] already
    # points at (SetFileLocations.py), so RunBootstrapsDCP.
    # BuildRealizationPayload's subsequent read needs zero changes.
    os.makedirs(GeneralDict['DCPJobDir'], exist_ok=True)
    with open(GeneralDict['FixtureFile'], 'w') as f:
        json.dump(result['fixtureJson'], f)

    return model


# ---------------------------------------------------------------------------
# 5. Top-level entry: called by FullSingleGalaxyFit.GalaxyFit() instead of
#    RW.RunWRKP(GeneralDict, GalaxyDict, 0) + RWF.ReadWRKPOutputFile when
#    GalaxyDict['UseDCP'] == 1.
# ---------------------------------------------------------------------------
def RunInitialFit(GeneralDict, GalaxyDict):
    cube_wcs = GetCubeWCS(GalaxyDict['CubeNameU'])
    payload = BuildInitialFitPayload(GeneralDict, GalaxyDict)

    result = LaunchInitialFitJob(GeneralDict, payload, objName=GalaxyDict.get('ObjName', ''))

    model = AdaptInitialFitResult(result, cube_wcs, GeneralDict, GalaxyDict)
    if not model.get('FITAchieved'):
        print("No best fit model made (JS initial fit)")
        if result.get('error'):
            print("  error:", result['error'])
    return model
