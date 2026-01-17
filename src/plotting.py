"""Bokeh visualization utilities for X-ray timing analysis."""

import numpy as np
from astropy.time import Time, TimeDelta
from bokeh.models import ColumnDataSource, Range1d, LinearAxis


def set_axis_styles(fig, font_size='14pt'):
    """Set consistent styling for Bokeh figures.

    Args:
        fig: Bokeh figure to style
        font_size: font size for labels (default '14pt')
    """
    fig.xaxis.axis_label_text_font_size = font_size
    fig.xaxis.major_label_text_font_size = font_size
    fig.yaxis.axis_label_text_font_size = font_size
    fig.yaxis.major_label_text_font_size = font_size

    fig.xaxis.major_label_text_font = 'Georgia'
    fig.xaxis.axis_label_text_font = 'Georgia'
    fig.yaxis.major_label_text_font = 'Georgia'
    fig.yaxis.axis_label_text_font = 'Georgia'

    fig.xaxis.axis_label_text_font_style = 'normal'
    fig.yaxis.axis_label_text_font_style = 'normal'

    if fig.title:
        fig.title.align = 'center'
        fig.title.text_font = 'Georgia'
        fig.title.text_font_size = font_size
        fig.title.text_font_style = 'bold'


def plot_lightcurve(p, lc, time_format='mjd', mjd_ref=None, gti=None, color='darkviolet', line_width=2):
    """Add lightcurve to Bokeh plot with gap handling.

    Args:
        p: Bokeh figure
        lc: Stingray Lightcurve object
        time_format: 'mjd', 'utc', or 'seconds'
        mjd_ref: reference MJD for 'seconds' format
        gti: GTI array (uses lc.gti if None)
        color: line color
        line_width: line width

    Returns: the figure
    """
    if time_format not in ['mjd', 'utc', 'seconds']:
        raise ValueError(f"time_format must be 'mjd', 'utc', or 'seconds', got '{time_format}'")

    if time_format == 'seconds' and mjd_ref is None:
        mjd_ref = lc.mjdref
        print(f"Using lightcurve MJD reference: {mjd_ref}")

    if gti is None and hasattr(lc, 'gti') and lc.gti is not None:
        gti = lc.gti

    lc_time = np.asarray(lc.time, dtype=np.float64)
    lc_countrate = np.asarray(lc.countrate, dtype=np.float64)

    def convert_time(time_sec):
        if time_format == 'mjd':
            return lc.mjdref + time_sec / 86400.0
        elif time_format == 'utc':
            t = Time(lc.mjdref, format='mjd') + TimeDelta(time_sec, format='sec')
            return t.to_datetime()
        else:
            return ((lc.mjdref + time_sec / 86400.0) - mjd_ref) * 86400.0

    if gti is not None:
        gti_array = np.asarray(gti, dtype=np.float64)
        for gti_start, gti_stop in gti_array:
            mask = (lc_time >= gti_start) & (lc_time <= gti_stop)
            if not np.any(mask):
                continue
            x_values = convert_time(lc_time[mask])
            source = ColumnDataSource(data=dict(x=x_values, y=lc_countrate[mask]))
            p.step('x', 'y', source=source, line_width=line_width, color=color, mode='center')
    else:
        x_values = convert_time(lc_time)
        y_values = np.where(lc_countrate == 0, np.nan, lc_countrate)
        source = ColumnDataSource(data=dict(x=x_values, y=y_values))
        p.step('x', 'y', source=source, line_width=line_width, color=color, mode='center')

    return p


def plot_folded_lightcurve(p, phase, countrate, err, n_bins, color='firebrick',
                           match_axis_color=False, scale_factor=None, line_alpha=1):
    """Plot folded lightcurve with error bars.

    Args:
        p: Bokeh figure
        phase: phase bin edges from fold() (n_bins+1)
        countrate: count rate per bin
        err: count rate errors
        n_bins: number of phase bins
        color: line and error bar color
        match_axis_color: match y-axis label to data color
        scale_factor: scale y-range based on max error
        line_alpha: line transparency
    """
    phase = phase[1:] - np.diff(phase) / 2

    extended_phase = np.concatenate([phase, phase + 1])
    extended_rate = np.concatenate([countrate, countrate])
    extended_err = np.concatenate([err, err])

    if scale_factor is not None:
        rate_min = np.min(extended_rate) - scale_factor * np.max(extended_err)
        rate_max = np.max(extended_rate) + scale_factor * np.max(extended_err)
        p.y_range.start = rate_min
        p.y_range.end = rate_max

    source = ColumnDataSource(data=dict(
        x=extended_phase,
        y=extended_rate,
        yerr=extended_err,
        err_bottom=extended_rate - extended_err,
        err_top=extended_rate + extended_err
    ))
    alpha_ratio = 0.4

    p.step('x', 'y', source=source, alpha=line_alpha, line_width=2, color=color, mode='center')
    p.vbar('x', width=1 / n_bins, top='err_top', bottom='err_bottom',
           source=source, color=color, line_color=None, alpha=line_alpha * alpha_ratio)

    if match_axis_color:
        p.yaxis.axis_label_text_color = color
        p.yaxis.major_label_text_color = color


def plot_hardness_ratio(p, phase, soft_countrate, soft_err, hard_countrate, hard_err, n_bins,
                        hr_color='black', scale_factor=6.0, line_width=2, line_alpha=1):
    """Plot hardness ratio on secondary y-axis.

    Args:
        p: Bokeh figure
        phase: phase bin edges (n_bins+1)
        soft_countrate, soft_err: soft band rates and errors
        hard_countrate, hard_err: hard band rates and errors
        n_bins: number of phase bins
        hr_color: HR line color
        scale_factor: HR axis range scaling
        line_width: line width
        line_alpha: line transparency

    Returns: (HR, HR_error) arrays
    """
    phase = phase[1:] - np.diff(phase) / 2

    HR = (hard_countrate - soft_countrate) / (hard_countrate + soft_countrate)
    HR_error = (2 * np.sqrt(hard_countrate**2 * soft_err**2 + soft_countrate**2 * hard_err**2)) / (hard_countrate + soft_countrate)**2

    extended_phase = np.concatenate([phase, phase + 1])
    extended_HR = np.concatenate([HR, HR])
    extended_HR_error = np.concatenate([HR_error, HR_error])

    hr_source = ColumnDataSource(data=dict(
        x=extended_phase,
        hr=extended_HR,
        hr_err_bottom=extended_HR - extended_HR_error,
        hr_err_top=extended_HR + extended_HR_error
    ))

    hr_min = np.min(HR) - scale_factor * 0.1
    hr_max = np.max(HR) + scale_factor * 0.1

    p.extra_y_ranges['hr'] = Range1d(hr_min, hr_max)

    alpha_ratio = 0.5
    p.step('x', 'hr', source=hr_source, line_width=2, alpha=line_alpha, color=hr_color,
           mode='center', y_range_name='hr')
    p.vbar('x', width=1 / n_bins, top='hr_err_top', bottom='hr_err_bottom',
           source=hr_source, color=hr_color, line_color=None, alpha=line_alpha * alpha_ratio, y_range_name='hr')

    ax2 = LinearAxis(y_range_name="hr", axis_label="Hardness Ratio")
    ax2.axis_label_text_color = hr_color
    ax2.major_label_text_color = hr_color
    ax2.axis_label_text_font_size = "12pt"
    ax2.axis_label_text_font_style = 'normal'
    ax2.major_label_text_font_size = "12pt"
    ax2.major_label_text_font = 'Georgia'
    ax2.axis_label_text_font = 'Georgia'
    p.add_layout(ax2, 'right')

    return HR, HR_error
