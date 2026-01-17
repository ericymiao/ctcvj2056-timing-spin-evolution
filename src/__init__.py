"""
CTCV J2056-3014 X-ray Timing Analysis Package

Modules:
    timing_analysis: Core periodogram and folding functions
    plotting: Bokeh visualization utilities
    statistics: Statistical tests for timing analysis
    demodulation: Orbital demodulation functions
"""

from .timing_analysis import (
    bayesian_periodogram,
    z2_confidence_interval,
    montgomery_frequency_uncertainty,
    get_frequency_uncertainty,
    fold,
    fold_lightcurve,
    phase_exposure_profile,
)

from .plotting import (
    set_axis_styles,
    plot_lightcurve,
    plot_folded_lightcurve,
    plot_hardness_ratio,
)

from .statistics import (
    calculate_hardness_ratio,
    test_hardness_ratio_variability,
    test_phase_correlation,
    print_observation_summary,
)

from .demodulation import (
    delta_t,
    shifted_times,
    demodulate_grid_search,
)

__version__ = '1.0.0'
__author__ = 'Eric Miao'
