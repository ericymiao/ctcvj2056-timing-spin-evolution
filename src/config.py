"""Paths and constants shared by the pipeline scripts and the notebooks.

Set the environment variable CTCV_DATA_DIR to the directory holding the cleaned, barycentered
event files (layout below). Everything else is relative to the repository.

    <CTCV_DATA_DIR>/
        analysis_files_strictfilter/   NICER: ni<obsid>_0mpu7_cl_bary.evt, the HENDRICS event
                                       files *_nicer_xti_ev_calib.nc and SCORPEON background
                                       light curves ni<obsid>mpu7_{0.3-2,2-emax}_bg.lc
        xmm/<obsid>/                   XMM: {pn,mos1,mos2}_clean_{src,bkg}.evt (barycentered)
        nustar/30801006002/event_cl/   NuSTAR: nu30801006002{A,B}01_cl_bary_{src,bkg}.evt
"""

import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_ROOT = os.environ.get('CTCV_DATA_DIR', '/Users/ericmiao/Documents/Kaya Research/CTCV_J2056')
NICER_DIR = f'{DATA_ROOT}/analysis_files_strictfilter'
XMM_DIR = f'{DATA_ROOT}/xmm'
NUSTAR_DIR = f'{DATA_ROOT}/nustar/30801006002/event_cl'

# Repository folders
DATA_DIR = os.path.join(REPO, 'data')
EVENTS_DIR = os.path.join(DATA_DIR, 'events')      # products of pipeline/01_prepare_events.py
RESULTS_DIR = os.path.join(REPO, 'results')        # numbers quoted in the paper (JSON)
FIGURE_DIR = os.path.join(REPO, 'paper_figures')   # figures of the paper

# ---------------------------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------------------------
NICER_OBS = {'4592010200': '2021A', '4592010300': '2021B', '7656010100': '2024'}
NICER_COMBINED = '4592010000'                       # 2021A + 2021B merged (combined periodogram)
XMM_OBS = {'0842570101': '2019', '0902500101': '2022'}

# Barycentered FITS event files used for the TOAs (PN files are SAS barycen products)
EVENT_FILES = {
    'nicer2021A': f'{NICER_DIR}/ni4592010200_0mpu7_cl_bary.evt',
    'nicer2021B': f'{NICER_DIR}/ni4592010300_0mpu7_cl_bary.evt',
    'nicer2024': f'{NICER_DIR}/ni7656010100_0mpu7_cl_bary.evt',
    'xmm2019_pn': f'{XMM_DIR}/0842570101/pn_clean_src_bary.evt',
    'xmm2022_pn': f'{XMM_DIR}/0902500101/pn_clean_src.evt',
}

# HENDRICS event files and SCORPEON background light curves read by the notebooks
NICER_NC = {obs: f'{NICER_DIR}/ni{obs}_0mpu7_cl_bary_nicer_xti_ev_calib.nc' for obs in list(NICER_OBS) + [NICER_COMBINED]}
NICER_BKG_LC = {obs: {'soft': f'{NICER_DIR}/ni{obs}mpu7_0.3-2_bg.lc', 'hard': f'{NICER_DIR}/ni{obs}mpu7_2-emax_bg.lc'} for obs in NICER_OBS}
NICER_EMAX = {'4592010200': 7, '4592010300': 10, '4592010000': 10, '7656010100': 10}   # keV; 2021A is cut at 7 keV
XMM_EVT = {obs: {module: {region: f'{XMM_DIR}/{obs}/{module}_clean_{region}.evt' for region in ['src', 'bkg']}
                 for module in ['pn', 'mos1', 'mos2']} for obs in XMM_OBS}
XMM_BKG_AREA_FACTOR = {'0842570101': (27**2) / (357**2), '0902500101': (20**2) / (113**2)}   # source / background extraction areas
NUSTAR_EVT = {fpm: {region: f'{NUSTAR_DIR}/nu30801006002{fpm[-1]}01_cl_bary_{region}.evt' for region in ['src', 'bkg']} for fpm in ['FPMA', 'FPMB']}
NUSTAR_BKG_AREA_FACTOR = (40**2) / (200**2)

# ---------------------------------------------------------------------------------------------
# Timing constants and pulsar-timing files
# ---------------------------------------------------------------------------------------------
P_SPIN = 29.60968584                 # spin period (s) used for folding and for the orbital search
F_SPIN = 1.0 / P_SPIN
FOLD_EPOCH_NICER = 59230             # MJD of phase zero for the folded NICER profiles
FOLD_EPOCH_XMM_NUSTAR = 58780
NICER_MJDREF = 56658 + 0.000777592592592593

TEMPLATE = os.path.join(DATA_DIR, 'template.gauss')                        # five-Gaussian template (nitemplate)
TEMPLATE_PHASES = os.path.join(DATA_DIR, 'nicer2021A_template_phases.npy')  # photon phases the template was fitted to
PAR_INITIAL = os.path.join(DATA_DIR, 'photon_toa_initial_model.par')        # model used by photon_toa
TOAS = os.path.join(DATA_DIR, 'toas.tim')                                   # one TOA per observation
TOAS_WITHIN = os.path.join(DATA_DIR, 'toas_within_2021A.tim')               # three TOAs within NICER 2021A
CRCUT_GTI = os.path.join(DATA_DIR, 'gti', 'nicer_2021A_crcut.ecsv')         # count-rate-cut GTI for the within-observation TOAs

# NICERsoft (git submodule) and the Python environment
NICERSOFT_SCRIPTS = os.path.join(REPO, 'external', 'NICERsoft', 'scripts')
