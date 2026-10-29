"""The figures of the paper (matplotlib). Each function returns the Figure; the notebooks save it.

Fonts and tick styling follow the paper: `paper_style()` returns the rcParams used by every
figure, and each function applies them through `plt.rc_context`. Pass `path=` to save the figure:
tick labels are created when the figure is drawn, so the file must be written inside the same
context for the styling to be complete.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter, FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker


def paper_style(fontsize=13):
    """rcParams used throughout the paper (inward ticks on all sides, constrained layout)."""
    return {'font.size': fontsize, 'axes.labelsize': fontsize, 'xtick.labelsize': fontsize, 'ytick.labelsize': fontsize,
            'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
            'xtick.minor.visible': True, 'ytick.minor.visible': True, 'figure.constrained_layout.use': True}


# ---------------------------------------------------------------------------------------------
# Periodograms (Figure 1 and Figure 8)
# ---------------------------------------------------------------------------------------------

def plot_periodogram(curves, figsize=(5, 3.6), xfmt=None, path=None):
    """Z^2 periodogram(s) in mHz.

    curves: list of dicts with 'frequency' (Hz), 'power', 'color', 'lw' and optionally
            'span' = (f_lo, f_hi, alpha) drawn as a shaded frequency-uncertainty band.
    xfmt:   e.g. '%.3f' to fix the tick label format.
    """
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=figsize)
        for c in curves:
            ax.plot(c['frequency'] * 1e3, c['power'], lw=c.get('lw', 1), color=c.get('color'))
            if c.get('span') is not None:
                f_lo, f_hi, alpha = c['span']
                ax.axvspan(f_lo * 1e3, f_hi * 1e3, alpha=alpha, color='darkviolet')
        ax.set_xlabel('Frequency (mHz)')
        ax.set_ylabel('Power')
        ax.ticklabel_format(axis='x', useOffset=False)
        if xfmt is not None:
            ax.xaxis.set_major_formatter(FormatStrFormatter(xfmt))
        if path is not None:
            fig.savefig(path)
    return fig


def plot_hires_periodogram(frequency, power, power_threshold, span_uhz=(33772.27, 33773.19),
                           xlim_uhz=(33771.75, 33773.75), ylim=(1050, 1800), path=None):
    """High-resolution combined periodogram in microhertz with the 2-sigma power threshold (Figure 8, top right)."""
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=(6.5, 4.3))
        ax.plot(frequency * 1e6, power, lw=1.5)
        ax.axvspan(span_uhz[0], span_uhz[1], alpha=0.15, color='darkviolet')
        ax.axhline(power_threshold, alpha=0.6, color='darkviolet', lw=2, ls='--')
        ax.set_xlim(*xlim_uhz)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r'Frequency ($\mu\text{Hz}$)')
        ax.set_ylabel('Power')
        ax.ticklabel_format(axis='x', useOffset=True)
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=True))
        if path is not None:
            fig.savefig(path)
    return fig


# ---------------------------------------------------------------------------------------------
# Folded profiles (Figures 2 and 3)
# ---------------------------------------------------------------------------------------------

def plot_folded_profiles(profiles, lw=2, alpha=0.35, ylim_top=None, figsize=(5, 3.6), path=None):
    """Step profiles over two cycles with shaded 1-sigma bands.

    profiles: list of (ext_phase, ext_rate, ext_err, color) already extended to [0, 2].
    """
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=figsize)
        for ext_phase, ext_rate, ext_err, color in profiles:
            ax.step(ext_phase, ext_rate, where='mid', color=color, lw=lw)
            ax.fill_between(ext_phase, ext_rate - ext_err, ext_rate + ext_err, step='mid', color=color, alpha=alpha)
        if ylim_top is not None:
            ax.set_ylim(top=ylim_top)
        ax.set_xlabel('Phase')
        ax.set_ylabel('Count Rate (cnts/s)')
        if path is not None:
            fig.savefig(path)
    return fig


def plot_folded_hr(ext_phase, total, total_err, hr, hr_err, path=None):
    """Background-subtracted total rate (red, left axis) with the hardness ratio (black, right axis) (Figure 3)."""
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=(5, 3.6))
        ext_total = np.concatenate([total, total])
        ext_total_err = np.concatenate([total_err, total_err])
        ax.step(ext_phase, ext_total, where='mid', color='firebrick', lw=2, alpha=0.5)
        ax.fill_between(ext_phase, ext_total - ext_total_err, ext_total + ext_total_err, step='mid', color='firebrick', alpha=0.08)
        ax.set_xlabel('Phase')
        ax.set_ylabel('Count Rate (cnts/s)', color='firebrick')
        ax.tick_params(axis='y', labelcolor='firebrick')
        ax.set_ylim(np.min(total) - 2 * np.max(total_err), np.max(total) + 2 * np.max(total_err))

        ext_hr = np.concatenate([hr, hr])
        ext_hr_err = np.concatenate([hr_err, hr_err])
        ax2 = ax.twinx()
        ax2.step(ext_phase, ext_hr, where='mid', color='black', lw=2)
        ax2.fill_between(ext_phase, ext_hr - ext_hr_err, ext_hr + ext_hr_err, step='mid', color='black', alpha=0.45)
        ax2.set_ylabel('Hardness Ratio')
        ax2.set_ylim(np.min(hr) - 2 * 0.1, np.max(hr) + 2 * 0.1)
        if path is not None:
            fig.savefig(path)
    return fig


# ---------------------------------------------------------------------------------------------
# Template (Figure 5)
# ---------------------------------------------------------------------------------------------

def plot_template(phases, template, nbins=50, path=None):
    """Folded photon phases (blue) with the Gaussian template (red) and fractional residuals (green)."""
    bins = np.linspace(0, 1, nbins + 1)
    bcs = (bins[1:] + bins[:-1]) / 2.0
    cod = np.histogram(phases, bins=bins, density=True)[0]
    err = (cod * float(nbins) / len(phases)) ** 0.5
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=(7 / 1.2, 5 / 1.2))

        # Fractional residuals on the right axis, drawn first so that they sit behind the profile
        ax_resid = ax.twinx()
        model = template(bcs)
        ax_resid.errorbar(bcs, cod / model - 1, yerr=err / model, color='green', ls=' ', marker='o', alpha=0.2)
        ax_resid.axhline(0, color='green', lw=1, alpha=0.7)
        ax_resid.set_ylabel(r'Fractional Residuals ($\times 10^{-2}$)', color='green')
        ax_resid.tick_params(axis='y', labelcolor='green')
        ax_resid.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x*100:.0f}'))
        ax_resid.set_xlim(0, 1)

        ax.step(bcs, cod, where='mid', color='blue', lw=1)
        ax.fill_between(bcs, cod - err, cod + err, step='mid', color='blue', alpha=0.15)
        dom = np.linspace(0, 1, 1000)
        ax.plot(dom, template(dom), color='red', lw=2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0.65, 1.4)
        ax.set_xlabel('Phase')
        ax.set_ylabel('Profile Amplitude')
        ax.set_zorder(ax_resid.get_zorder() + 1)
        ax.patch.set_visible(False)
        if path is not None:
            fig.savefig(path)
    return fig


# ---------------------------------------------------------------------------------------------
# Pulsar timing residuals (Figures 6, 7, 8 bottom, 10, 11) and the initial-frequency chart (Figure 9)
# ---------------------------------------------------------------------------------------------

def plot_toa_residuals(mjds, resids, errs, zero_line=False, ylim=None, path=None):
    """Residuals of the five TOAs, one panel (Figure 10)."""
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=(7, 2.5))
        ax.scatter(mjds, resids, marker='x', s=150, linewidths=2.5, color='firebrick', zorder=5)
        ax.bar(mjds, height=2 * errs, bottom=resids - errs, width=60, alpha=0.2, color='firebrick', edgecolor='none')
        if zero_line:
            ax.axhline(0, color='gray', ls='--')
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.set_xlabel('MJD')
        ax.set_ylabel('Residuals (s)')
        if path is not None:
            fig.savefig(path)
    return fig


def plot_prefit_montage(panels, path=None):
    """Three stacked pre-fit residual panels for successive initial models (Figure 6).

    panels: list of (label, mjds, resids, errs).
    """
    with plt.rc_context(paper_style(12)):
        fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True)
        for ax, (label, t, r, e) in zip(axes, panels):
            ax.scatter(t, r, marker='x', s=150, linewidths=2.5, color='firebrick', zorder=5)
            ax.bar(t, height=2 * e, bottom=r - e, width=60, alpha=0.2, color='firebrick', edgecolor='none')
            span = np.max(r) - np.min(r)
            ax.set_ylim(np.min(r) - 0.15 * span, np.max(r) + 0.35 * span)
            ax.set_ylabel('Residuals (s)')
            ax.text(0.02, 0.92, label, transform=ax.transAxes, fontsize=11, va='top')
        axes[-1].set_xlabel('MJD')
        if path is not None:
            fig.savefig(path)
    return fig


def plot_within_residuals(mjds, resids, errs, bar_width, ylim=(-0.2, 0.2), path=None):
    """Residuals of the three within-observation TOAs, one small panel (Figure 7)."""
    with plt.rc_context(paper_style(12)):
        fig, ax = plt.subplots(figsize=(4 / 1.2, 3 / 1.2))
        ax.set_ylim(*ylim)
        ax.scatter(mjds, resids, marker='x', s=200, linewidths=3, zorder=5)
        ax.bar(mjds, height=2 * errs, bottom=resids - errs, width=bar_width, alpha=0.2, edgecolor='none', color='#1f77b4')
        ax.axhline(0, color='gray', ls='--', lw=0.8)
        ax.set_xlabel('MJD')
        ax.set_ylabel('Residuals (s)')
        if path is not None:
            fig.savefig(path)
    return fig


def plot_resid_slopes_montage(panels, bar_width, path=None):
    """Within-observation pre-fit residuals for the five candidate aliases (Figure 8, bottom).

    panels: dict label -> (mjds, resids, errs, f0) for 'leftpeak2', 'leftpeak', 'rightpeak', 'rightpeak2', 'center'.
    """
    with plt.rc_context(paper_style(12)):
        fig = plt.figure(figsize=(9, 5.2))
        gs = GridSpec(2, 3, figure=fig)
        positions = {'leftpeak2': gs[0, 0], 'leftpeak': gs[0, 1], 'rightpeak': gs[1, 0], 'rightpeak2': gs[1, 1], 'center': gs[:, 2]}
        for label, pos in positions.items():
            mjds, resids, errs, f0 = panels[label]
            ax = fig.add_subplot(pos)
            ax.scatter(mjds, resids, marker='x', s=200, linewidths=3, zorder=5)
            ax.bar(mjds, height=2 * errs, bottom=resids - errs, width=bar_width, alpha=0.2, edgecolor='none', color='#1f77b4')
            ax.axhline(0, color='gray', ls='--', lw=0.8)
            ax.set_ylim(-1.5, 1.5)
            ax.set_xlabel('MJD')
            ax.set_ylabel('Residuals (s)')
            ax.set_title(f'{label}: $f$ = {f0:.8f} Hz', fontsize=11)
        if path is not None:
            fig.savefig(path)
    return fig


def plot_error_ellipse(corr, sigmas, x_scale=1e-10, y_scale=1e-18, path=None):
    """1, 2 and 3 sigma joint confidence ellipses of (F0, F1) with the marginal 1-sigma bands (Figure 11)."""
    from scipy import linalg, stats
    vals, vecs = linalg.eigh(corr)
    with plt.rc_context(paper_style()):
        fig, ax = plt.subplots(figsize=(5.7, 4.3))
        colors = {1: 'salmon', 2: 'indianred', 3: 'maroon'}
        for n_sigma in [3, 2, 1]:
            thresh = np.sqrt(stats.chi2(2).isf(2 * stats.norm.cdf(-n_sigma)))
            angles = np.linspace(0, 2 * np.pi, 200)
            points = thresh * (np.sqrt(vals[0]) * np.cos(angles)[:, None] * vecs[None, :, 0]
                               + np.sqrt(vals[1]) * np.sin(angles)[:, None] * vecs[None, :, 1])
            ax.plot(points[:, 0] * sigmas[0] / x_scale, points[:, 1] * sigmas[1] / y_scale,
                    label=f'{n_sigma} sigma', color=colors[n_sigma], lw=2)
        ax.axvspan(-sigmas[0] / x_scale, sigmas[0] / x_scale, alpha=0.2, color='cornflowerblue')
        ax.axhspan(-sigmas[1] / y_scale, sigmas[1] / y_scale, alpha=0.2, color='cornflowerblue')
        ax.set_xlabel(r'$\Delta$ F0 ($\times 10^{-10}$ Hz)')
        ax.set_ylabel(r'$\Delta$ F1 ($\times 10^{-18}$ Hz/s)')
        ax.legend(loc='upper right')
        if path is not None:
            fig.savefig(path)
    return fig


def plot_initial_frequency(labels, freqs, uncerts, delta_f, path=None):
    """Successive initial-frequency estimates with their uncertainties and the tolerable error band (Figure 9)."""
    F_final = freqs[-1]
    with plt.rc_context(dict(matplotlib.rcParamsDefault)):
        fig, ax = plt.subplots(figsize=(13, 3))
        ax.axvspan(F_final - delta_f, F_final + delta_f, color='C2', alpha=0.15, label=r'Max tolerable $\Delta f$')
        ax.axvline(F_final, color='k', linestyle='--', linewidth=1.5)
        y_pos = np.arange(len(labels))[::-1]
        for i in range(len(labels)):
            ax.barh(y_pos[i], width=2 * uncerts[i], left=freqs[i] - uncerts[i], height=0.4, color='C0', alpha=0.3)
        ax.scatter(freqs, y_pos, marker='x', color='C0', s=200, zorder=5, linewidths=2.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=13)
        ax.set_xlabel(r'$F$ (Hz)', fontsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.grid(axis='x', alpha=0.3)
        fig.tight_layout()
        if path is not None:
            fig.savefig(path, bbox_inches='tight', dpi=300)
    return fig
