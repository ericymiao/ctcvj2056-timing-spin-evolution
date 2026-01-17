"""Core timing analysis: Z^2 periodograms and phase folding."""

import numpy as np
from scipy.stats import ncx2
from scipy.special import erf
from astropy.time import Time, TimeDelta
from stingray.gti import cross_gtis


# ========================= Z^2 PERIODOGRAM FUNCTIONS =========================

def bayesian_periodogram(frequencies, powers, nharm):
    """Convert Z^2 powers to Bayesian posterior probabilities.

    Args:
        frequencies: frequency array
        powers: Z^2_n power array
        nharm: number of harmonics

    Returns: posterior probability array
    """
    noise_variance = 2 * nharm
    return np.exp(powers / noise_variance)


def z2_confidence_interval(peak_power, nharm, sigma_level):
    """Get power confidence interval using noncentral chi-squared.

    Args:
        peak_power: detected Z^2_n power
        nharm: number of harmonics
        sigma_level: confidence level (e.g. 1.0 for 1-sigma)

    Returns: (power_lower, power_upper)
    """
    dof = 2 * nharm
    nc_param = peak_power - dof
    central_prob = erf(sigma_level / np.sqrt(2))
    lower_tail = (1 - central_prob) / 2
    upper_tail = 1 - lower_tail

    power_lower = ncx2.ppf(lower_tail, dof, nc_param)
    power_upper = ncx2.ppf(upper_tail, dof, nc_param)
    return power_lower, power_upper


def _get_boundaries_from_level(x, y, level, x0):
    """Find where y crosses level around x0."""
    max_idx = np.argmin(np.abs(x - x0))
    idx = max_idx
    min_x = max_x = x0

    while idx > 0 and y[idx] > level:
        min_x = x[idx]
        idx -= 1

    idx = max_idx
    while idx < y.size and y[idx] > level:
        max_x = x[idx]
        idx += 1

    return min_x, max_x


def montgomery_frequency_uncertainty(peak_freq, peak_power, nharm, obs_time, n_events, sigma_level=1.0):
    """Frequency uncertainty using Montgomery method.

    Uses formula: sigma_f = (sqrt(6) / (pi*T)) / sqrt(2*(z_max - 2m))
    Assumes strong detection (Z^2 >> 2n).

    Args:
        peak_freq: peak frequency (Hz)
        peak_power: peak Z^2_n power
        nharm: number of harmonics
        obs_time: total observation time (seconds)
        n_events: number of events
        sigma_level: confidence level (default 1.0)

    Returns: (freq_lower, freq_upper, freq_uncertainty)
    """
    if peak_power <= 2 * nharm:
        raise ValueError(f"Peak power ({peak_power:.1f}) must be >> 2*nharm ({2*nharm})")
    if obs_time <= 0 or n_events <= 0:
        raise ValueError("obs_time and n_events must be positive")

    delta_chi2 = peak_power - 2 * nharm
    sigma_freq = sigma_level * (np.sqrt(6.0) / (np.pi * obs_time)) * (1.0 / np.sqrt(2.0 * delta_chi2))
    freq_lower, freq_upper = peak_freq - sigma_freq, peak_freq + sigma_freq

    peak_period = 1.0 / peak_freq
    period_lower, period_upper = 1.0 / freq_upper, 1.0 / freq_lower
    period_err_upper, period_err_lower = period_upper - peak_period, peak_period - period_lower

    print("=" * 70)
    print(f"MONTGOMERY METHOD ({sigma_level}-sigma):")
    print(f"  Detection: Z^2_n = {peak_power:.1f}, delta_chi2 = {delta_chi2:.1f}")
    print(f"  T = {obs_time:.1f} s, m = {nharm}")
    print(f"  Frequency: {peak_freq:.15e} +/- {sigma_freq:.6e} Hz")
    print(f"  Period:    {peak_period:.15e} +{period_err_upper:.6e} -{period_err_lower:.6e} s")
    print("=" * 70)

    return freq_lower, freq_upper, sigma_freq


def get_frequency_uncertainty(frequencies, powers, peak_freq, peak_power, nharm, sigma_level):
    """Get frequency uncertainty by tracing periodogram at confidence level.

    Args:
        frequencies: periodogram frequency array
        powers: periodogram power array
        peak_freq: peak frequency
        peak_power: peak power value
        nharm: number of harmonics
        sigma_level: confidence level in sigma

    Returns: (freq_min, freq_max, power_lower, power_upper)
    """
    power_lower, power_upper = z2_confidence_interval(peak_power, nharm, sigma_level)
    freq_min, freq_max = _get_boundaries_from_level(frequencies, powers, power_lower, peak_freq)

    peak_freq = float(peak_freq)
    freq_min = float(freq_min)
    freq_max = float(freq_max)
    peak_power = float(peak_power)
    power_lower = float(power_lower)
    power_upper = float(power_upper)

    freq_err_upper = freq_max - peak_freq
    freq_err_lower = peak_freq - freq_min

    peak_period = 1 / peak_freq
    period_min = 1 / freq_max
    period_max = 1 / freq_min
    period_err_lower = peak_period - period_min
    period_err_upper = period_max - peak_period

    print("=" * 70)
    print("FREQUENCY RESULTS:")
    print(f"  Peak:        {peak_freq:.15e} Hz")
    print(f"  Uncertainty: +{freq_err_upper:.6e} / -{freq_err_lower:.6e} Hz")
    print()
    print("PERIOD RESULTS:")
    print(f"  Peak:        {peak_period:.15e} s")
    print(f"  Uncertainty: +{period_err_upper:.6e} / -{period_err_lower:.6e} s")
    print("=" * 70)

    return freq_min, freq_max, power_lower, power_upper


# ========================= FOLDING FUNCTIONS =========================

def _to_seconds_since_epoch(times, mjdref, mjd_epoch, time_is_mjd=True):
    """Convert times to seconds since epoch."""
    time_mjd = times if time_is_mjd else mjdref + times / 86400.0
    return (time_mjd - mjd_epoch) * 86400.0


def _to_phases(time_seconds, period):
    """Convert seconds to phase (0 to 1)."""
    return (time_seconds / period) % 1.0


def _compute_gti_exposure(gtis_sec, period, first_cycle, last_cycle, n_bins):
    """Compute exposure per phase bin from GTIs."""
    exposure = np.zeros(n_bins)
    phase_edges = np.linspace(0, 1, n_bins + 1)

    for i in range(n_bins):
        cycles = np.arange(first_cycle - 1, last_cycle + 1)
        phase_gtis = np.column_stack([
            (cycles + phase_edges[i]) * period,
            (cycles + phase_edges[i + 1]) * period
        ])
        exposure[i] = np.sum(np.diff(cross_gtis((phase_gtis, gtis_sec))))

    return exposure


def _distribute_contributions(start, end, rate, period, n_bins, counts, exposure):
    """Distribute lightcurve bin across phase bins (in-place)."""
    phase_start, phase_end = (start / period) % 1, (end / period) % 1
    n_periods = int(end / period) - int(start / period)
    phase_bins = np.linspace(0, 1, n_bins + 1)

    def add_contribution(p_start, p_end):
        start_idx = np.clip(int(p_start * n_bins), 0, n_bins - 1)
        end_idx = np.clip(int(p_end * n_bins), 0, n_bins - 1)

        if start_idx == end_idx:
            overlap = (p_end - p_start) * period
            exposure[start_idx] += overlap
            counts[start_idx] += rate * overlap
        else:
            for i in range(start_idx, min(end_idx + 1, n_bins)):
                overlap_start = max(p_start, phase_bins[i])
                overlap_end = min(p_end, phase_bins[i + 1] if i < n_bins - 1 else 1.0)
                if overlap_end > overlap_start:
                    overlap = (overlap_end - overlap_start) * period
                    exposure[i] += overlap
                    counts[i] += rate * overlap

    if n_periods == 0:
        add_contribution(phase_start, phase_end)
    else:
        add_contribution(phase_start, 1.0)
        if n_periods > 1:
            time_per_bin = (n_periods - 1) * period / n_bins
            counts[:] += rate * time_per_bin
            exposure[:] += time_per_bin
        add_contribution(0.0, phase_end)


def phase_exposure_profile(events, period, mjd_epoch, n_bins):
    """Calculate exposure per phase bin from event GTIs.

    Args:
        events: Stingray EventList with gti, mjdref
        period: folding period (seconds)
        mjd_epoch: MJD of phase zero
        n_bins: number of phase bins

    Returns: exposure array (seconds per bin)
    """
    event_gti = np.asarray(events.gti, dtype=np.float64)
    gtis = Time(events.mjdref, format='mjd') + TimeDelta(event_gti, format='sec') - Time(mjd_epoch, format='mjd')
    gtis_sec = gtis.to_value('sec')

    first_cycle = int(gtis_sec[0, 0] // period)
    last_cycle = int(gtis_sec[-1, 1] // period)

    return _compute_gti_exposure(gtis_sec, period, first_cycle, last_cycle, n_bins)


def fold(events, period, mjd_epoch, n_bins):
    """Fold events into phase bins.

    Args:
        events: Stingray EventList with time, gti, mjdref
        period: folding period (seconds)
        mjd_epoch: MJD of phase zero
        n_bins: number of phase bins

    Returns: (bins, counts, countrate, err, exposure)
    """
    event_time = np.asarray(events.time, dtype=np.float64)
    time_sec = _to_seconds_since_epoch(event_time, events.mjdref, mjd_epoch, time_is_mjd=False)
    phases = _to_phases(time_sec, period)

    bins = np.linspace(0, 1, n_bins + 1)
    counts, _ = np.histogram(phases, bins)
    exposure = phase_exposure_profile(events, period, mjd_epoch, n_bins)

    epsilon = 1e-10
    countrate = counts / (exposure + epsilon)
    err = np.sqrt(counts) / (exposure + epsilon)

    return bins, counts, countrate, err, exposure


def fold_lightcurve(lightcurve, period, mjd_epoch, n_bins):
    """Fold binned lightcurve with exposure correction.

    Args:
        lightcurve: Stingray Lightcurve with time, countrate, mjdref, dt
        period: folding period (seconds)
        mjd_epoch: MJD of phase zero
        n_bins: number of phase bins

    Returns: (bins, counts, countrate, err, exposure)
    """
    times = np.asarray(lightcurve.time, dtype=np.float64)
    countrates = np.asarray(lightcurve.countrate, dtype=np.float64)
    mjdref = lightcurve.mjdref
    bin_width = lightcurve.dt

    has_counts = hasattr(lightcurve, 'counts') and lightcurve.counts is not None
    if has_counts:
        lc_counts = np.asarray(lightcurve.counts, dtype=np.float64)
        has_counts = np.sum(lc_counts) > 0

    if has_counts:
        lc_exposure = np.where(countrates > 0, lc_counts / countrates, bin_width)
        lc_exposure = np.where((countrates == 0) & (lc_counts > 0), bin_width, lc_exposure)
    else:
        lc_counts = countrates * bin_width
        lc_exposure = np.full_like(times, bin_width, dtype=np.float64)

    time_sec = _to_seconds_since_epoch(times, mjdref, mjd_epoch, time_is_mjd=False)
    bin_starts = time_sec - bin_width / 2
    bin_ends = time_sec + bin_width / 2

    phase_counts = np.zeros(n_bins)
    phase_exposure = np.zeros(n_bins)

    n_periods = (bin_ends // period).astype(int) - (bin_starts // period).astype(int)
    simple_mask = n_periods == 0

    if np.any(simple_mask):
        simple_phases = _to_phases(bin_starts[simple_mask], period)
        simple_indices = np.clip((simple_phases * n_bins).astype(int), 0, n_bins - 1)
        np.add.at(phase_exposure, simple_indices, lc_exposure[simple_mask])
        np.add.at(phase_counts, simple_indices, lc_counts[simple_mask])

    for idx in np.where(~simple_mask)[0]:
        rate = countrates[idx]
        _distribute_contributions(bin_starts[idx], bin_ends[idx], rate,
                                 period, n_bins, phase_counts, phase_exposure)

    bins = np.linspace(0, 1, n_bins + 1)
    epsilon = 1e-10
    countrate = phase_counts / (phase_exposure + epsilon)
    err = np.sqrt(phase_counts) / (phase_exposure + epsilon)

    return bins, phase_counts, countrate, err, phase_exposure
