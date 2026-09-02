"""Orbital demodulation functions for X-ray timing analysis."""

import numpy as np
from multiprocessing import Pool, Value, Lock
from stingray.pulse import z_n_search

# Module-level shared state for parallel progress tracking
_progress_counter = None
_progress_lock = None
_progress_total = None
_progress_step = None


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


def make_T_grid(T_min, T_max, T_tot):
    """Generate T grid with spacing DeltaT = T^2/T_tot (uniform in 1/T).

    This ensures independence between adjacent sampling points in
    Fourier wave number space.

    Args:
        T_min: minimum period (seconds)
        T_max: maximum period (seconds)
        T_tot: total observation time span (seconds)

    Returns: 1D array of T values from T_min to T_max
    """
    u = np.arange(1/T_min, 1/T_max, -1/T_tot)
    return 1.0 / u


def _init_progress(counter, lock, total, step):
    """Initialize shared progress state in each worker."""
    global _progress_counter, _progress_lock, _progress_total, _progress_step
    _progress_counter = counter
    _progress_lock = lock
    _progress_total = total
    _progress_step = step


def _worker_grid_search(args):
    """Worker function for parallel grid search over a chunk of T values."""
    T_chunk, event_times, P_spin, As, phis, nharm = args
    results = np.zeros((len(T_chunk), len(As), len(phis)))
    for i_T, T in enumerate(T_chunk):
        for i_A, A in enumerate(As):
            for i_phi, phi in enumerate(phis):
                t_shifted = shifted_times(phi, A, event_times, T)
                _, z = z_n_search(t_shifted, 1 / P_spin, nharm=nharm)
                results[i_T, i_A, i_phi] = z
        # Update shared progress counter
        with _progress_lock:
            _progress_counter.value += 1
            done = _progress_counter.value
            if done % _progress_step == 0 or done == _progress_total:
                pct = 100 * done / _progress_total
                print(f"{pct:.0f}% ({done}/{_progress_total} T values)", flush=True)
    return results


def demodulate_grid_search_parallel(event_times, P_spin, Ts, As, phis,
                                    nharm=1, n_workers=8):
    """Parallelized grid search over orbital parameters to maximize Z^2 power.

    Splits the T array across workers. Each worker evaluates all (A, phi)
    combinations for its assigned T values.

    Args:
        event_times: event arrival times (seconds)
        P_spin: spin period (seconds)
        Ts: orbital periods to search (seconds)
        As: timing delay amplitudes to search (seconds)
        phis: orbital phases to search (radians)
        nharm: harmonics for Z^2 test
        n_workers: number of parallel workers

    Returns: 3D array of Z^2 powers, shape (len(Ts), len(As), len(phis))
    """
    t = np.asarray(event_times)
    total = len(Ts)
    step = max(1, total // 50)  # print every ~2%
    counter = Value('i', 0)
    lock = Lock()
    T_chunks = np.array_split(Ts, n_workers)
    args = [(chunk, t, P_spin, As, phis, nharm) for chunk in T_chunks]

    with Pool(n_workers, initializer=_init_progress,
              initargs=(counter, lock, total, step)) as pool:
        chunk_results = pool.map(_worker_grid_search, args)

    return np.concatenate(chunk_results, axis=0)
