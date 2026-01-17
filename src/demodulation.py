"""Orbital demodulation functions for X-ray timing analysis."""

import numpy as np
from stingray.pulse import z_n_search


def delta_t(phi, A, t, T):
    """Calculate orbital timing delay.

    Args:
        phi: orbital phase offset (radians)
        A: timing delay amplitude (seconds)
        t: event times (seconds)
        T: orbital period (seconds)

    Returns: time delay array
    """
    return A * np.sin(2 * np.pi * t / T - phi)


def shifted_times(phi, A, times, T):
    """Apply orbital correction to event times.

    Args:
        phi: orbital phase offset (radians)
        A: timing delay amplitude (seconds)
        times: event times (seconds)
        T: orbital period (seconds)

    Returns: corrected times with orbital delay removed
    """
    deltas = delta_t(phi, A, times, T)
    return times - deltas

def shifted_times_scramble(phi, A, times, T):
    """Randomly scramble correction times to events, for significance simulations. 

    Args:
        phi: orbital phase offset (radians)
        A: timing delay amplitudes (seconds)
        times: event times (seconds)
        T: orbital period (seconds)
    """
    deltas = np.random.shuffle(delta_t(phi, A, times, T).copy())
    return times - deltas


def demodulate_grid_search(event_times, P_spin, Ts, As, phis, nharm=2, verbose=True):
    """Grid search over orbital parameters to maximize Z^2 power.

    Args:
        event_times: event arrival times (seconds)
        P_spin: spin period (seconds)
        Ts: orbital periods to search (seconds)
        As: timing delay amplitudes to search (seconds)
        phis: orbital phases to search (radians)
        nharm: harmonics for Z^2 test
        verbose: print progress

    Returns: 3D array of Z^2 powers, shape (len(Ts), len(As), len(phis))
    """
    t = np.asarray(event_times)
    results = np.zeros((len(Ts), len(As), len(phis)))

    total_T = len(Ts)

    for i_T, T in enumerate(Ts):
        for i_A, A in enumerate(As):
            for i_phi, phi in enumerate(phis):
                t_shifted = shifted_times(phi, A, t, T)
                _, z = z_n_search(t_shifted, 1 / P_spin, nharm=nharm)
                results[i_T, i_A, i_phi] = z

        if verbose:
            print(f"{100 * (i_T + 1) / total_T:.1f}%")

    return results


def find_best_parameters(results, Ts, As, phis):
    """Find orbital parameters that maximize Z^2 power.

    Args:
        results: 3D array from demodulate_grid_search
        Ts, As, phis: searched parameter arrays

    Returns: dict with best T, A, phi and z2_max
    """
    max_idx = np.unravel_index(np.argmax(results), results.shape)
    i_T, i_A, i_phi = max_idx

    return {
        'T': Ts[i_T],
        'A': As[i_A],
        'phi': phis[i_phi],
        'z2_max': results[max_idx],
        'indices': max_idx
    }


def project_to_period_axis(results, Ts):
    """Project demodulation results onto orbital period axis.

    For each period, finds max Z^2 over all amplitude/phase combinations.

    Args:
        results: 3D array from demodulate_grid_search
        Ts: orbital periods

    Returns: max Z^2 power at each period
    """
    max_powers = np.zeros(len(Ts))
    for i in range(len(Ts)):
        max_powers[i] = np.max(results[i, :, :])
    return max_powers
