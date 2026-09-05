"""Event-list loaders and observation summaries (the numbers of Table 1)."""

import numpy as np
from stingray import Lightcurve
import hendrics.io as HENio

import config


# ---------------------------------------------------------------------------------------------
# NICER
# ---------------------------------------------------------------------------------------------

def load_nicer_events(obs):
    """Return {'total', 'soft', 'hard'} EventLists for a NICER ObsID (0.3-emax, 0.3-2, 2-emax keV)."""
    emax = config.NICER_EMAX[obs]
    total = HENio.load_events(config.NICER_NC[obs])
    total.filter_energy_range([0.3, emax], inplace=True)
    return {'total': total, 'soft': total.filter_energy_range([0.3, 2]), 'hard': total.filter_energy_range([2, emax])}


def load_nicer_background(obs):
    """SCORPEON background light curves {'soft', 'hard'} for a NICER ObsID."""
    return {band: Lightcurve.read(config.NICER_BKG_LC[obs][band], rate_column='countrate', err_column='BACK_ERR')
            for band in ['soft', 'hard']}


def nicer_summary(events, bkg):
    """Good time, span, counts, mean rates and background-subtracted (net) counts of one observation.

    The background counts are the SCORPEON model rates (0.5 s bins) summed over the observation.
    """
    ev = events['total']
    good_time = float(np.sum(ev.gti[:, 1] - ev.gti[:, 0]))
    span = float(ev.time[-1] - ev.time[0])
    total_counts = int(len(ev.time))
    bkg_counts = float(np.sum(bkg['soft'].countrate * 0.5) + np.sum(bkg['hard'].countrate * 0.5))
    lc = ev.to_lc(dt=20)
    lc.apply_gtis(inplace=True)
    return {'good_time_ks': good_time / 1000, 'span_ks': span / 1000, 'total_counts': total_counts,
            'mean_rate_cps': float(np.mean(lc.countrate)),
            'mean_bkg_rate_cps': float(np.mean(bkg['soft'].countrate) + np.mean(bkg['hard'].countrate)),
            'bkg_counts': bkg_counts, 'net_counts': total_counts - bkg_counts}


# ---------------------------------------------------------------------------------------------
# XMM-Newton
# ---------------------------------------------------------------------------------------------

def load_xmm_events(obs):
    """Nested dict events[module][region][band] for one XMM ObsID (0.3-10, 0.3-2, 2-10 keV)."""
    events = {}
    for module in ['pn', 'mos1', 'mos2']:
        events[module] = {}
        for region in ['src', 'bkg']:
            total = HENio.load_events(config.XMM_EVT[obs][module][region])
            total.filter_energy_range([0.3, 10], inplace=True)
            events[module][region] = {'total': total, 'soft': total.filter_energy_range([0.3, 2]),
                                      'hard': total.filter_energy_range([2, 10])}
    return events


def xmm_summary(events, obs):
    """Per-module exposure, source, area-scaled background and net counts (Table 1)."""
    factor = config.XMM_BKG_AREA_FACTOR[obs]
    out = {}
    for module in events:
        src = len(events[module]['src']['total'].time)
        bkg = len(events[module]['bkg']['total'].time) * factor
        gti = events[module]['src']['total'].gti
        out[module] = {'good_time_ks': float(np.sum(gti[:, 1] - gti[:, 0])) / 1000, 'source_counts': int(src),
                       'bkg_counts_scaled': float(bkg), 'net_counts': float(src - bkg)}
    return out


# ---------------------------------------------------------------------------------------------
# NuSTAR
# ---------------------------------------------------------------------------------------------

def load_nustar_events():
    """Nested dict events[fpm][region][band] (3-79, 3-10, 10-79 keV)."""
    events = {}
    for fpm in ['FPMA', 'FPMB']:
        events[fpm] = {}
        for region in ['src', 'bkg']:
            total = HENio.load_events(config.NUSTAR_EVT[fpm][region])
            total.filter_energy_range([3, 79], inplace=True)
            events[fpm][region] = {'total': total, 'soft': total.filter_energy_range([3, 10]),
                                   'hard': total.filter_energy_range([10, 79])}
    return events


def nustar_summary(events):
    """Per-module exposure, source, area-scaled background and net counts (Table 1)."""
    out = {}
    for fpm in events:
        src = len(events[fpm]['src']['total'].time)
        bkg = len(events[fpm]['bkg']['total'].time) * config.NUSTAR_BKG_AREA_FACTOR
        gti = events[fpm]['src']['total'].gti
        out[fpm] = {'good_time_ks': float(np.sum(gti[:, 1] - gti[:, 0])) / 1000, 'source_counts': int(src),
                    'bkg_counts_scaled': float(bkg), 'net_counts': float(src - bkg)}
    return out


# ---------------------------------------------------------------------------------------------
# Folded profiles
# ---------------------------------------------------------------------------------------------

def phase_centers(phase_edges):
    """Bin centres of a folded profile, repeated over two cycles for plotting."""
    centers = phase_edges[1:] - np.diff(phase_edges) / 2
    return centers, np.concatenate([centers, centers + 1])


def subtract_background(profile, bkg_profile, factor=1.0):
    """Background-subtract a folded profile (phase, counts, rate, err, exposure) tuple."""
    phase, counts, rate, err, exposure = profile
    _, _, bkg_rate, bkg_err, _ = bkg_profile
    return phase, counts, rate - bkg_rate * factor, np.sqrt(err**2 + (bkg_err * factor)**2), exposure

