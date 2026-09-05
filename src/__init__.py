"""
CTCV J2056-3014 X-ray timing analysis package

Modules:
    config:          data paths, observation table and timing constants
    data_io:         event-list loaders and observation summaries (Table 1)
    timing_analysis: Z^2 periodograms, peak uncertainties and phase folding
    statistics:      hardness-ratio tests
    demodulation:    orbital Roemer-delay grid search
    plotting:        the figures of the paper
    events:          FITS event-file utilities for the TOA pipeline
"""

__version__ = '2.0.0'
__author__ = 'Eric Miao'
