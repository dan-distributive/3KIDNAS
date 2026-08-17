"""
RunBootstrapsDCP.py
===================

Single-dispatch DCP bootstrap driver: one DCP job, one slice per bootstrap
realization, each slice running the ENTIRE per-realization pipeline native
Fortran runs atomically (Bootstrap_Error_Analysis.GetBootstrapModel):
resample -> SoFiA -> InitialAnalysis -> fit. Mirrors
js/app/bootstrap-realization-launcher.js.

Replaces the prior two-round-trip design (bootstrap-resample-launcher.js +
local RunFixtureOnlyFit + bootstrap-fit-launcher.js): that design existed
only because the ~2500-line Fortran InitialAnalysis chain (moment maps,
shape/radial-profile estimation, noise) had no JS port and had to run
natively, locally, between two separate DCP dispatches. Now that it's
ported (js/src/PreAnalysis/), that local step -- and the
second DCP dispatch -- are both gone.

NOT YET RUN END-TO-END: the JS side (bootstrap-realization-launcher.js and
everything it calls) has been extensively unit- and integration-tested
locally. This Python file has NOT been exercised against a real DCP
dispatch or a real Fortran-produced fixture -- review carefully, and start
with a small nBootstraps run_both.js-style trial before trusting it at
production scale.
"""

import os
import json
import copy
import subprocess
import time
import base64
import warnings
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.wcs import FITSFixedWarning

from . import Bootstrap_Outputs as BO

# Fortran PipelineGlobals.f hardcoded parameter constants -- not read from
# any config file, so there's nothing to parse; just replicate the literals.
V_ROT_LIMS = [0.0, 400.0]
SIZE_LIMS = [2.0, 10.0]
LIN_SD_LIMS_MSOL = [0.1, 50.0]


# ---------------------------------------------------------------------------
# 1. Resampling geometry (unchanged from the prior two-job design).
# ---------------------------------------------------------------------------
def ComputeBsCent(GalaxyDict):
    """
    Resampling-geometry centre/PA/inc, in the units genFlipBootstrapSample
    (the JS port) expects. Mirrors MakeBootstrapSample.WriteBootstrapFile's
    own GeoStr calculation EXACTLY (same formula, duplicated rather than
    shared, so this file doesn't touch that still-used, already-proven
    function -- keep the two in sync if that formula ever changes).
    """
    Model = GalaxyDict['BestFitModel']
    RefVel = GalaxyDict['CubeHeader']['CRVAL3']
    RefChan = GalaxyDict['CubeHeader']['CRPIX3']
    dV = GalaxyDict['CubeHeader']['CDELT3'] / 1000.
    DeltaV = Model['VSYS'][0] - RefVel / 1000.
    VCenter = DeltaV / dV + RefChan
    return {
        'centX': float(Model['XCENTER'][0]),
        'centY': float(Model['YCENTER'][0]),
        'centV': float(VCenter),
        'pa': float((Model['POSITIONANGLE'][0] + 90.) * np.pi / 180.),
        'inc': float(Model['INCLINATION'][0] * np.pi / 180.),
    }


# ---------------------------------------------------------------------------
# 2. Config values the fitting-options file carries that DumpFittingFixture
#    does NOT dump (checked directly against SingleGalaxyFitTests.f's
#    DumpFittingFixture -- it stops at nRingsPerBeam). Extracted by
#    replaying FittingOptionsInputs.f's own sequential read order over the
#    SAME options-file lines Python already loads (GeneralDict
#    ['MainWRKPOptionsLines'], via RunWRKP.LoadDefaultWRKPFiles) -- no
#    Fortran change needed.
#
#    nTargRings/RadGrid are read here too (to walk past them correctly) but
#    NOT used -- for a bootstrap realization the ring grid is always pinned
#    to the ORIGINAL fit's own best-fit radii (RunWRKP.WriteWRKPOptionsFile's
#    BootstrapSwitch==1 branch: nRTarg=len(BestFitModel['R']),
#    RArr=BestFitModel['R'], computed at dispatch time, not read from the
#    generic template), so BuildRealizationPayload builds nTargRings/
#    radGridArcsec from GalaxyDict['BestFitModel']['R'] instead.
# ---------------------------------------------------------------------------
def ParseFittingOptionsExtras(optionsLines):
    it = iter(optionsLines)

    def read_value():
        next(it)          # label/comment line
        return next(it).strip()

    core_code_switch = int(read_value())        # noqa: F841 (walked, not used)
    likelihood_switch = int(read_value())        # noqa: F841
    param_interface_switch = int(read_value())   # noqa: F841
    moment_map_switch = int(read_value())

    if moment_map_switch == 3:
        next(it); next(it); next(it)             # label, SecondaryCubeFile, SecondaryMaskFile
    elif moment_map_switch == 2:
        next(it); next(it); next(it); next(it); next(it)  # label, 3 MapNames, ProfileName

    center_source = int(read_value())
    shape_source = int(read_value())             # noqa: F841
    size_source = int(read_value())              # noqa: F841
    vsys_source = int(read_value())              # noqa: F841
    linear_log_sd_switch = int(read_value())     # noqa: F841 (also in the fixture; walked for position only)
    cmode = int(read_value())                    # noqa: F841 (also in the fixture)
    cloud_base_surf_dens = float(read_value())   # noqa: F841 (also in the fixture)
    sigma_lengths = float(read_value())          # noqa: F841 (also in the fixture, as observedBeam.sigmaLengths)
    noise_sigma_lim = float(read_value())
    n_targ_rings_raw = int(read_value())         # noqa: F841 -- see docstring: not used, bootstrap grid is pinned instead

    if n_targ_rings_raw > 0:
        next(it)  # label
        next(it)  # RadStr line (not used here)

    return {
        'centerSource': center_source,
        'noiseSigmaLim': noise_sigma_lim,
    }


def ComputeSDLims(pixelSizeX, pixelSizeY, channelSize):
    """
    Mirrors PreAnalysis.f's PreGalaxyAnalysis preamble (SD limits block,
    lines ~26-39): derives LinSDLims/LogSDLims from the hardcoded
    LinSDLims_MSol (PipelineGlobals.f) and this cube's own pixel/channel
    geometry. Not read from any config file -- these are DERIVED, so
    there's nothing to parse, just the same formula replicated.
    """
    JyAS_To_MsolPC = 1.24756e20 / (6.0574e5 * 1.823e18 * (2.0 * np.pi / np.log(256.0)))
    sd_conv1 = 1.0 / abs(pixelSizeX * pixelSizeY)
    lin_sd_lims = []
    log_sd_lims = []
    for v_msol in LIN_SD_LIMS_MSOL:
        sd_jy_as2 = v_msol * JyAS_To_MsolPC       # MSolPc2_To_JyAS2
        sd = sd_jy_as2 / sd_conv1
        sd = sd / abs(channelSize)
        lin_sd_lims.append(sd)
        log_sd_lims.append(float(np.log10(sd)))
    return lin_sd_lims, log_sd_lims


# ---------------------------------------------------------------------------
# 3. Build the single shared payload every realization's slice uses
#    identically (resample geometry + fit config). Per-realization variation
#    comes only from realizationIndex itself (the resample step's own RNG
#    derivation), not from anything in this payload.
# ---------------------------------------------------------------------------
def BuildRealizationPayload(GeneralDict, GalaxyDict, fixture_path):
    with open(fixture_path, 'r') as f:
        fixture = json.load(f)

    bs_cent = ComputeBsCent(GalaxyDict)

    extras = ParseFittingOptionsExtras(GeneralDict['MainWRKPOptionsLines'])

    lin_sd_lims, log_sd_lims = ComputeSDLims(
        fixture['observedDC']['pixelSizeX'],
        fixture['observedDC']['pixelSizeY'],
        fixture['observedDC']['channelSize'],
    )

    # Pinned ring grid: bootstrap realizations always use the ORIGINAL fit's
    # own best-fit radii, not whatever the generic fitting-options file says
    # (matches RunWRKP.WriteWRKPOptionsFile's BootstrapSwitch==1 branch,
    # computed at dispatch time from BestFitModel, not read from a file).
    best_fit_radii = list(GalaxyDict['BestFitModel']['R'])

    # RAW (un-brightness-converted, Jy/beam) FITS bytes for resample+SoFiA.
    # fixture['observedDC']'s own validFlux is Jy/pixel -- dumped by
    # Fortran AFTER DCBrightnessConversion (SingleGalaxyFitTests.f:45-51,91),
    # which runs right after the cube is loaded, before PreGalaxyAnalysis/
    # SoFiA-catalogue consumption. But SoFiA ITSELF never sees that converted
    # cube: it's invoked as an external binary directly on GalaxyDict['CubeName']
    # (the observed cube) and GalaxyDict['BestFitModel']['ModelCube'] (the
    # best-fit model cube) via MakeBootstrapSample -> BootStrapSampler, whose
    # own LoadCubesForBootstrap has its DCBrightnessConversion calls
    # commented out (BootstrapRuntimeInputs.f:105-106) -- confirmed empirically
    # too: BootStrapSampler's own resampled-cube output sums to ~28x (its
    # own beamAreaPixels) the fixture's converted validFlux for the same
    # geometry/seed. Found via a direct resampled-cube-sum trace comparison
    # (native's own catalog kin_pa 66.454915 vs JS's previous
    # pre-converted-input catalog kin_pa 66.556271 -- a real ~0.1 deg PA
    # divergence, not stochastic spread, since Nelder-Mead is deterministic
    # and idum-in matched exactly at eval 1). Reading these here and
    # resampling/SoFiA-ing on the RAW bytes (JS applies the /beamAreaPixels
    # conversion itself, AFTER SoFiA, right before InitialAnalysis --
    # see bootstrap-realization-launcher.js) matches native's real order of
    # operations exactly instead of approximately.
    with open(GalaxyDict['CubeName'], 'rb') as f:
        observed_raw_fits_b64 = base64.b64encode(f.read()).decode('ascii')
    with open(GalaxyDict['BestFitModel']['ModelCube'], 'rb') as f:
        model_raw_fits_b64 = base64.b64encode(f.read()).decode('ascii')

    # observedDC header only -- NOT the dense validFlux/flattendValidIndices
    # arrays the fixture also carries under that key (the fixture's own
    # already-brightness-converted flux, unused now that resample+SoFiA runs
    # on the raw FITS bytes below instead -- see that comment). Confirmed
    # unread anywhere in runBootstrapRealization (it only pulls scalar header
    # fields off observedDC); stripping them, and dropping modelFlux entirely
    # (same story -- kept for a hypothetical future consumer, never read),
    # cuts the dispatched payload from ~9.4MB to ~1.4MB per job.
    observed_dc_header = {k: v for k, v in fixture['observedDC'].items()
                          if k not in ('validFlux', 'flattendValidIndices')}

    return {
        # ---- Resample+SoFiA geometry (same fields the old round-1 payload carried) ----
        'observedDC': observed_dc_header,
        'observedCubeRawFitsB64': observed_raw_fits_b64,  # RAW Jy/beam FITS bytes, base64 -- resample+SoFiA input
        'modelCubeRawFitsB64': model_raw_fits_b64,        # RAW Jy/beam FITS bytes, base64 -- resample+SoFiA input
        'beamMajorAxis': fixture['observedBeam']['beamMajorAxis'],
        'beamMinorAxis': fixture['observedBeam']['beamMinorAxis'],
        'bsCent': bs_cent,
        'velBlockSize': 1,
        'bootstrapSeed': GeneralDict.get('BootstrapSeed', 0),
        'sofiaParTemplate': ''.join(GeneralDict['SoFiATemplate']),

        # ---- Fit-side config (from the fixture + newly-parsed options-file
        #      extras + hardcoded PipelineGlobals.f constants) ----
        'observedBeam': fixture['observedBeam'],  # full beam (kernel/sigma vector) for the fit's own convolution
        'centerSource': extras['centerSource'],
        'nRingsPerBeam': fixture['nRingsPerBeam'],
        'nTargRings': len(best_fit_radii),
        'radGridArcsec': best_fit_radii,
        'sdSwitch': fixture['linearLogSDSwitch'],
        # Which of the 13 ring properties are constant-across-rings /
        # permanently fixed -- required by initialAnalysis() (see its own
        # comment on why there's no safe default). Found missing via a
        # direct numerical cross-check against Fortran's real ProbeSwitch
        # chi2: without these, the JS side silently fits 65 free params
        # (13 per ring x 5 rings, "everything free") instead of Fortran's
        # real 16, producing a garbage chi2 and (very likely) explaining
        # both the large Inc/PA/VSys divergences seen in the real DCP
        # bootstrap comparison and a meaningful share of the eval-count/
        # performance gap (a 66-vertex simplex instead of a 17-vertex one).
        'constParams': fixture['constParams'],
        'fixedParams': fixture['fixedParams'],
        'linSDLims': lin_sd_lims,
        'logSDLims': log_sd_lims,
        'vRotLims': V_ROT_LIMS,
        'sizeLims': SIZE_LIMS,
        'noiseSigmaLim': extras['noiseSigmaLim'],
        'cmode': fixture['cmode'],
        'cloudBaseSurfDens': fixture['cloudBaseSurfDens'],
        # Static per-galaxy fit seed -- Fortran's real per-realization
        # bootstrap fit never varies idum across realizations (RunWRKP /
        # WriteWRKPMainFile always reuses the template's static value), so
        # this is intentionally the SAME value for every realization, not
        # derived from realizationIndex.
        'fitIdum': fixture['idum'],
        'ftol': fixture['ftol'],
        'likelihoodSwitch': fixture['likelihoodSwitch'],
    }


def GetCubeWCS(cube_path):
    hdr = fits.getheader(cube_path)
    # Silences astropy's benign "datfix ... Set MJD-OBS" notice -- WCS's own
    # header-normalization fix, not a problem with the cube or the fit.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        w = wcs.WCS(hdr)
    return w


# ---------------------------------------------------------------------------
# 4. Launch the single DCP job -- one slice per realization.
# ---------------------------------------------------------------------------
def LaunchDCPRealizationJob(GeneralDict, payload, nBoot, objName=''):
    workdir = GeneralDict['DCPJobDir']
    launcher = GeneralDict['DCPRealizationLauncher']
    os.makedirs(workdir, exist_ok=True)

    payload_path = os.path.join(workdir, 'realization_payload.json')
    results_path = os.path.join(workdir, 'dcp_realization_results.json')

    if os.path.exists(results_path):
        os.remove(results_path)

    with open(payload_path, 'w') as f:
        json.dump(payload, f)

    cmd = [
        'node', launcher,
        '--nBootstraps', str(nBoot),
        '--payload', payload_path,
        '--out', results_path,
        '--jobName', '3KIDNAS-' + objName,
    ]
    # DCP_API_KEY unset -> fall back to the launcher's own --local N mode:
    # direct, in-process execution of runBootstrapRealization (no dcp-client,
    # no network, no credentials needed -- see the launcher's own --local
    # branch for why it deliberately avoids job.localExec()/dcp-client
    # entirely rather than using it for "local" execution). Real DCP
    # dispatch (--apiKey) is used whenever the credential is present.
    # DCP_FORCE_LOCAL=1: opt-in override for local/manual testing, so a run
    # with a real DCP_API_KEY already set in the environment (e.g. for some
    # other, real DCP dispatch elsewhere in the same session) can still
    # force THIS invocation to run locally instead of dispatching for real.
    # (Formerly also used by run_all_three.js's since-removed "js-local"
    # leg, which compared this local-execution path against real js-dcp
    # dispatch -- consistently ~4x slower, so that comparison leg is gone,
    # but the underlying local-execution override is still generally useful.)
    apiKey = os.environ.get('DCP_API_KEY')
    forceLocal = os.environ.get('DCP_FORCE_LOCAL') == '1'
    if apiKey and not forceLocal:
        cmd += ['--apiKey', apiKey]
    else:
        if not apiKey:
            print("DCP_API_KEY not set -- falling back to --local %d (direct in-process execution, no dispatch)" % nBoot)
        cmd += ['--local', str(nBoot)]
    computeGroups = os.environ.get('DCP_COMPUTE_GROUPS')
    if computeGroups:
        cmd += ['--computeGroups', computeGroups]
    slicePrice = os.environ.get('DCP_SLICE_PRICE')
    if slicePrice:
        cmd += ['--slicePrice', slicePrice]

    r = subprocess.run(cmd, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError("DCP bootstrap-realization job failed (returncode %d)" % r.returncode)

    with open(results_path, 'r') as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# 5. Adapter: worker returns -> BootstrapModels list-of-dicts. UNCHANGED from
#    the prior two-job design -- the merged launcher returns the exact same
#    per-realization result shape bootstrap-fit-launcher.js used to.
# ---------------------------------------------------------------------------
def AdaptResults(worker_results, cube_wcs, nRealizations):
    by_idx = {}
    for res in worker_results:
        by_idx[res['realizationIndex']] = res

    models = [None] * nRealizations
    for step in range(nRealizations):
        res = by_idx.get(step)
        if res is None or res.get('sofiaFailed') or not res.get('converged', False):
            models[step] = {'FITAchieved': False}
            continue

        xc = np.array(res['XCENTER'], dtype=float)
        yc = np.array(res['YCENTER'], dtype=float)
        try:
            ra, dec = cube_wcs.celestial.all_pix2world(xc, yc, 0)
        except Exception:
            ra = np.full_like(xc, np.nan)
            dec = np.full_like(yc, np.nan)

        model = {
            'FITAchieved': True,
            'XCENTER': list(xc),
            'YCENTER': list(yc),
            'INCLINATION': np.array(res['INCLINATION'], dtype=float),
            'POSITIONANGLE': np.array(res['POSITIONANGLE'], dtype=float),
            'VSYS': np.array(res['VSYS'], dtype=float),
            'VROT': np.array(res['VROT'], dtype=float),
            'VDISP': np.array(res['VDISP'], dtype=float),
            'SURFDENS': np.array(res['SURFDENS'], dtype=float),
            'SURFDENS_FACEON': np.array(res['SURFDENS'], dtype=float),
            'RA': list(ra),
            'DEC': list(dec),
            'R': np.array(res.get('R', list(range(len(res['VROT'])))), dtype=float),
            'R_SD': np.array(res.get('R_SD', list(range(len(res['SURFDENS'])))), dtype=float),
        }
        models[step] = model
    return models


def SummarizeWorkerTimings(results, round_label, dispatch_wall_time):
    exec_times = [r['timings']['totalMs'] for r in results if r.get('timings')]
    if not exec_times:
        return None
    sum_exec_s = sum(exec_times) / 1000.
    print("%s: worker execution min=%dms max=%dms mean=%dms sum=%.1fs"
          % (round_label, min(exec_times), max(exec_times),
             sum(exec_times) / len(exec_times), sum_exec_s))
    print("%s: DCP overhead (dispatch wall time - sum of worker execution time): %.1fs"
          % (round_label, dispatch_wall_time - sum_exec_s))
    return sum_exec_s


# ---------------------------------------------------------------------------
# 6. Top-level entry: called by FullSingleGalaxyFit.GalaxyFit() when
#    GalaxyDict['UseDCP'] == 1.
# ---------------------------------------------------------------------------
def RunBootstrapsDCP(GeneralDict, GalaxyDict):
    nBoot = GalaxyDict['nBootstraps']

    cube_wcs = GetCubeWCS(GalaxyDict['CubeName'])
    payload = BuildRealizationPayload(GeneralDict, GalaxyDict, GeneralDict['FixtureFile'])

    DispatchStart = time.time()
    worker_results = LaunchDCPRealizationJob(GeneralDict, payload, nBoot, objName=GalaxyDict.get('ObjName', ''))
    GalaxyDict['DCP_DispatchTime'] = time.time() - DispatchStart
    GalaxyDict['DCP_DispatchWorkerTimeSum'] = SummarizeWorkerTimings(
        worker_results, 'Bootstrap realizations (resample+SoFiA+fit)', GalaxyDict['DCP_DispatchTime'])

    # Kept for FullSingleGalaxyFit.py's existing timing-CSV logic, which
    # still expects a "resample" checkpoint distinct from "dispatch" --
    # with everything in one dispatch now, there's no separate wall-clock
    # phase to report here, so this is just 0 (informational only).
    GalaxyDict['DCP_ResampleTime'] = 0.0
    GalaxyDict['DCP_LocalPrepTime'] = 0.0

    BootstrapModels = AdaptResults(worker_results, cube_wcs, nBoot)

    BootstrapTimings = []
    worker_by_idx = {r['realizationIndex']: r for r in worker_results}
    for step in range(nBoot):
        r = worker_by_idx.get(step, {})
        t = r.get('timings') or {}
        BootstrapTimings.append({
            'realizationIndex': step,
            'resampleMs': t.get('resampleMs'),
            'sofiaMs': t.get('sofiaMs'),
            # Now measures the in-slice InitialAnalysis chain's own time,
            # not a separate local native subprocess call -- schema name
            # kept for run_both.js compatibility.
            'fixtureFitMs': t.get('fixtureFitMs'),
            'fitMs': t.get('fitMs'),
            'convolveMs': t.get('convolveMs'),
            'convolveCalls': t.get('convolveCalls'),
            'evalCount': t.get('evalCount'),
            'totalMs': t.get('totalMs'),
        })
    BO.StoreBootstrapTimings_JSON(GalaxyDict, BootstrapTimings)

    nSucc = sum(1 for m in BootstrapModels if m.get('FITAchieved'))
    print("DCP bootstrap: %d/%d realizations succeeded" % (nSucc, nBoot))
    print("DCP timing: dispatch %.1fs" % (GalaxyDict['DCP_DispatchTime'],))
    return BootstrapModels
