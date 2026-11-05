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
