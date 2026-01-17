"""Statistical tests for X-ray timing analysis."""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import f as f_dist
from scipy.stats import norm


def calculate_hardness_ratio(soft, hard, soft_err, hard_err):
    """Calculate HR = (H-S)/(H+S) with error propagation.

    Args:
        soft, hard: soft/hard band count rates
        soft_err, hard_err: rate errors

    Returns: (hr, hr_err) arrays
    """
    soft = np.asarray(soft)
    hard = np.asarray(hard)
    soft_err = np.asarray(soft_err)
    hard_err = np.asarray(hard_err)

    hr = (hard - soft) / (hard + soft)
    hr_err = (2 * np.sqrt(hard**2 * soft_err**2 + soft**2 * hard_err**2)) / (hard + soft)**2

    return hr, hr_err


def _sinusoid(phase, offset, amplitude, phase_shift):
    """Sinusoidal model for fitting."""
    return offset + amplitude * np.sin(2 * np.pi * phase + phase_shift)


def test_hardness_ratio_variability(hr_data, hr_errors, phase=None, verbose=True):
    """F-test comparing constant vs sinusoidal model for HR.

    Args:
        hr_data: hardness ratio values per phase bin
        hr_errors: errors on HR values
        phase: phase values 0-1 (uniform if None)
        verbose: print results

    Returns: dict with f_statistic, p_value, n_sigma, chi2 values, fit params
    """
    hr_data = np.asarray(hr_data)
    hr_errors = np.asarray(hr_errors)
    n_data = len(hr_data)

    if phase is None:
        phase = np.linspace(0, 1, n_data, endpoint=False) + 0.5 / n_data
    else:
        phase = np.asarray(phase)

    # Model 1: Constant (weighted mean)
    weights = 1.0 / hr_errors**2
    mean_hr = np.sum(weights * hr_data) / np.sum(weights)
    chi2_constant = np.sum(((hr_data - mean_hr) / hr_errors)**2)
    dof_constant = n_data - 1

    # Model 2: Sinusoid
    try:
        popt, pcov = curve_fit(_sinusoid, phase, hr_data,
                               p0=[mean_hr, np.std(hr_data), 0.0],
                               sigma=hr_errors, absolute_sigma=True)
        offset_fit, amplitude_fit, phase_shift_fit = popt
        perr = np.sqrt(np.diag(pcov))

        hr_model = _sinusoid(phase, *popt)
        chi2_sinusoid = np.sum(((hr_data - hr_model) / hr_errors)**2)
        dof_sinusoid = n_data - 3
    except RuntimeError:
        return {
            'f_statistic': np.nan,
            'p_value': 1.0,
            'n_sigma': 0.0,
            'result': 'Sinusoid fit failed'
        }

    # F-test
    delta_chi2 = chi2_constant - chi2_sinusoid
    delta_dof = dof_constant - dof_sinusoid
    F_statistic = (delta_chi2 / delta_dof) / (chi2_sinusoid / dof_sinusoid)
    p_value = 1.0 - f_dist.cdf(F_statistic, delta_dof, dof_sinusoid)

    # Convert to sigma
    if p_value > 0 and p_value < 1:
        n_sigma = norm.ppf(1 - p_value / 2)
    elif p_value <= 0:
        n_sigma = np.inf
    else:
        n_sigma = 0.0

    if p_value < 0.001:
        interpretation = "Sinusoid significantly better (p < 0.001)"
    elif p_value < 0.05:
        interpretation = "Sinusoid better (p < 0.05)"
    else:
        interpretation = "No significant improvement"

    if verbose:
        print(f"\nModel 1 (Constant): chi2 = {chi2_constant:.2f}, dof = {dof_constant}, chi2_red = {chi2_constant/dof_constant:.2f}")
        print(f"Model 2 (Sinusoid): chi2 = {chi2_sinusoid:.2f}, dof = {dof_sinusoid}, chi2_red = {chi2_sinusoid/dof_sinusoid:.2f}")
        print(f"  Offset: {offset_fit:.4f} +/- {perr[0]:.4f}")
        print(f"  Amplitude: {amplitude_fit:.4f} +/- {perr[1]:.4f}")
        print(f"  Phase shift: {phase_shift_fit/(2*np.pi):.3f} cycles")
        print(f"\nF-test: F = {F_statistic:.2f}, p = {p_value:.2e}, {n_sigma:.1f} sigma")
        print(f"Result: {interpretation}")

    return {
        'f_statistic': F_statistic,
        'p_value': p_value,
        'n_sigma': n_sigma,
        'chi2_constant': chi2_constant,
        'chi2_sinusoid': chi2_sinusoid,
        'fit_params': popt,
        'fit_errors': perr,
        'result': interpretation
    }


def test_phase_correlation(flux_data, flux_errors, hr_data, hr_errors, phase=None, verbose=True):
    """Test phase relationship between flux and hardness ratio.

    Fits sinusoids to both, compares phase of maximum.
    Anticorrelation: phase diff ~ 0.5; Correlation: phase diff ~ 0

    Args:
        flux_data, flux_errors: flux values and errors per bin
        hr_data, hr_errors: HR values and errors per bin
        phase: phase values 0-1 (uniform if None)
        verbose: print results

    Returns: dict with phase_diff, phase_max values, amplitude significances
    """
    flux_data = np.asarray(flux_data)
    flux_errors = np.asarray(flux_errors)
    hr_data = np.asarray(hr_data)
    hr_errors = np.asarray(hr_errors)
    n_bins = len(hr_data)

    if phase is None:
        phase = np.linspace(0, 1, n_bins, endpoint=False) + 0.5 / n_bins
    else:
        phase = np.asarray(phase)

    # Fit sinusoids
    try:
        popt_hr, pcov_hr = curve_fit(_sinusoid, phase, hr_data,
                                      p0=[np.mean(hr_data), np.std(hr_data), 0.0],
                                      sigma=hr_errors, absolute_sigma=True)
        offset_hr, amp_hr, phase_shift_hr = popt_hr
        err_hr = np.sqrt(np.diag(pcov_hr))

        popt_flux, pcov_flux = curve_fit(_sinusoid, phase, flux_data,
                                          p0=[np.mean(flux_data), np.std(flux_data), 0.0],
                                          sigma=flux_errors, absolute_sigma=True)
        offset_flux, amp_flux, phase_shift_flux = popt_flux
        err_flux = np.sqrt(np.diag(pcov_flux))
    except RuntimeError:
        return {
            'phase_diff': np.nan,
            'result': 'Sinusoid fit failed'
        }

    # Amplitude significance
    amp_sig_hr = abs(amp_hr) / err_hr[1]
    amp_sig_flux = abs(amp_flux) / err_flux[1]

    # Normalize amplitudes to positive
    if amp_hr < 0:
        amp_hr = -amp_hr
        phase_shift_hr += np.pi
    if amp_flux < 0:
        amp_flux = -amp_flux
        phase_shift_flux += np.pi

    # Phase of maximum
    phase_max_hr = ((np.pi / 2 - phase_shift_hr) / (2 * np.pi)) % 1
    phase_max_flux = ((np.pi / 2 - phase_shift_flux) / (2 * np.pi)) % 1
    err_phase_max_hr = err_hr[2] / (2 * np.pi)
    err_phase_max_flux = err_flux[2] / (2 * np.pi)

    # Phase difference
    delta_phase = phase_max_hr - phase_max_flux
    if delta_phase > 0.5:
        delta_phase -= 1.0
    elif delta_phase < -0.5:
        delta_phase += 1.0

    err_delta_phase = np.sqrt(err_phase_max_hr**2 + err_phase_max_flux**2)

    # Significance tests
    deviation_from_anticorr = abs(abs(delta_phase) - 0.5)
    significance_anticorr = deviation_from_anticorr / err_delta_phase
    significance_from_zero = abs(delta_phase) / err_delta_phase

    if significance_anticorr < significance_from_zero and significance_anticorr < 3.0:
        interpretation = "ANTICORRELATED - HR max aligns with flux min"
    elif significance_from_zero < 3.0:
        interpretation = "IN PHASE - HR max aligns with flux max"
    else:
        interpretation = "Intermediate phase shift"

    if verbose:
        print(f"\nHR amplitude: {amp_hr:.4f} +/- {err_hr[1]:.4f} ({amp_sig_hr:.1f} sigma)")
        print(f"Flux amplitude: {amp_flux:.4f} +/- {err_flux[1]:.4f} ({amp_sig_flux:.1f} sigma)")
        print(f"\nFlux max at phase: {phase_max_flux:.3f} +/- {err_phase_max_flux:.4f}")
        print(f"HR max at phase: {phase_max_hr:.3f} +/- {err_phase_max_hr:.4f}")
        print(f"Phase difference: {delta_phase:.3f} +/- {err_delta_phase:.4f} cycles")
        print(f"\nDeviation from in-phase (delta_phi=0): {significance_from_zero:.1f} sigma")
        print(f"Deviation from anticorr (|delta_phi|=0.5): {significance_anticorr:.1f} sigma")
        print(f"Result: {interpretation}")

    return {
        'phase_diff': delta_phase,
        'phase_diff_err': err_delta_phase,
        'phase_max_flux': phase_max_flux,
        'phase_max_hr': phase_max_hr,
        'amp_significance_flux': amp_sig_flux,
        'amp_significance_hr': amp_sig_hr,
        'result': interpretation
    }


def print_observation_summary(events, bkg_lightcurves=None, obs_list=None, instrument='NICER'):
    """Print formatted observation summary.

    Args:
        events: dict of EventList objects keyed by obs ID
        bkg_lightcurves: dict of background Lightcurves (optional)
        obs_list: list of obs IDs to summarize (default: all)
        instrument: instrument name for header
    """
    if obs_list is None:
        obs_list = list(events.keys())

    print(f"{instrument} OBSERVATION SUMMARY")

    for obs in obs_list:
        if obs not in events:
            continue

        evt = events[obs]
        if hasattr(evt, 'total'):
            evt = evt['total']

        total_time = np.sum(evt.gti[:, 1] - evt.gti[:, 0])
        total_counts = len(evt.time)

        print(f"\n OBSERVATION: {obs}")
        print("-" * 50)
        print(f"Total Time:     {total_time/1000:8.2f} ks")
        print(f"Total Counts:    {total_counts:8,}")

        if bkg_lightcurves is not None and obs in bkg_lightcurves:
            bkg = bkg_lightcurves[obs]
            if 'soft' in bkg and 'hard' in bkg:
                total_background = (np.sum(bkg['soft'].countrate * 0.5) +
                                   np.sum(bkg['hard'].countrate * 0.5))
                avg_background = (np.mean(bkg['soft'].countrate) +
                                 np.mean(bkg['hard'].countrate))
                source_counts = total_counts - total_background

                print(f"Avg Bkg Rate:    {avg_background:8.3f} cts/s")
                print(f"Total Bkg:       {total_background:8.0f}")
                print(f"Source Counts:   {source_counts:8,.0f}")

    print()
