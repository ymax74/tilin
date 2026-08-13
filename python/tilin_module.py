#!/usr/bin/env python3
"""
Pandas-based API for TILIN132 visual-binary orbit fitting.

Wraps the algorithmic core in tilin.py (Thiele-Innes method, ported from the
Pulkovo Observatory Pascal program TILIN132) so it can be called from other
Python scripts, using pandas DataFrames instead of input.txt / output.txt.

Usage:
    import pandas as pd
    from tilin_module import fit_orbit_dataframe

    obs = pd.DataFrame({
        'year':  [1836.21, 1852.92, 1857.90, 1880.59, 1891.84],
        'rho':   [2.50,    2.89,    2.63,    2.10,    2.04],
        'theta': [295.60,  292.50,  292.10,  295.70,  293.00],
    })

    result = fit_orbit_dataframe(obs, year_col='year', rho_col='rho', theta_col='theta')
    print(result.elements)   # {'a':..., 'e':..., 'T':..., 'P':..., 'omega':..., 'Node':..., 'i':...}
    print(result.residuals)  # DataFrame of O-C residuals, one row per input observation

    # Bootstrap uncertainty estimate + corner plot:
    samples = bootstrap_orbit_elements(result, n_boot=100)
    limits  = orbit_element_limits(samples)
    plot_orbit_corner(samples, best_fit=result.elements, filename='orbit_corner.png')
    plot_orbit_family(samples, result, filename='orbit_family.png')
"""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

import tilin as _core

ORBIT_ELEMENTS = ('a', 'e', 'P', 'T', 'omega', 'Node', 'i')

ELEMENT_LABELS = {
    'a':     'a (arcsec)',
    'e':     'e',
    'T':     'T (yr)',
    'P':     'P (yr)',
    'omega': 'ω (deg)',
    'Node':  'Ω (deg)',
    'i':     'i (deg)',
}


@dataclass
class OrbitFitResult:
    elements: dict          # Campbell elements: a (arcsec), e, T (yr), P (yr), omega/Node/i (deg)
    t0: float                # reference epoch (yr), the observations' mean epoch
    rho0: float               # separation at t0 (arcsec)
    theta0: float             # position angle at t0 (deg)
    rho_dot: float            # dρ/dt at t0 (arcsec/yr)
    theta_dot: float          # dθ/dt at t0 (deg/yr)
    curvature: float          # radius of curvature of the apparent path at t0 (arcsec)
    rv: float                 # line-of-sight velocity at t0 (arcsec/yr, parallax=1)
    mass_pi1: float           # total mass a³/P² (solar masses, parallax=1)
    rms: float                # weighted RMS residual over the used (non-excluded) observations
    rms_all: float            # weighted RMS residual over all observations
    n_obs: int                # number of input observations
    n_excluded: int           # number of observations rejected as outliers
    residuals: pd.DataFrame   # per-observation O-C table, indexed like the input DataFrame


def fit_orbit_dataframe(df, year_col='year', rho_col='rho', theta_col='theta', verbose=True):
    """
    Fit a visual-binary orbit (Thiele-Innes / TILIN132 method) from a
    DataFrame of relative astrometric observations.

    df        — DataFrame with (at least) the three columns below. Rows with
                NaN in any of the three are dropped before fitting.
    year_col  — column with the epoch of observation, in years.
    rho_col   — column with the angular separation, in arcsec.
    theta_col — column with the position angle, in degrees.
    verbose   — print fit-progress messages (global search / outlier
                rejection rounds) to stdout.

    Returns an OrbitFitResult. Outliers are detected and excluded the same
    way the original program does (iterative sigma-clipping); flagged rows
    are kept in .residuals with excluded=True rather than dropped, so the
    caller can inspect or override the decision.
    """
    work = df[[year_col, rho_col, theta_col]].dropna()
    if len(work) < 3:
        raise ValueError("need at least 3 observations (with non-NaN year/rho/theta) to fit an orbit")

    epochs_abs = work[year_col].astype(float).tolist()
    rho_obs    = work[rho_col].astype(float).tolist()
    theta_obs  = work[theta_col].astype(float).tolist()

    ti, excluded, t0, epochs_rel = _core.fit_orbit(epochs_abs, rho_obs, theta_obs, verbose=verbose)

    weights = [1.0] * len(epochs_abs)
    summary = _core.compute_fit_summary(ti, excluded, t0, epochs_rel, epochs_abs,
                                        rho_obs, theta_obs, weights)

    residuals = pd.DataFrame(summary['rows'], index=work.index)

    elements = {k: summary[k] for k in ('a', 'e', 'T', 'P', 'omega', 'Node', 'i')}

    return OrbitFitResult(
        elements=elements,
        t0=summary['t0'],
        rho0=summary['rho0'],
        theta0=summary['theta0'],
        rho_dot=summary['rho_dot'],
        theta_dot=summary['theta_dot'],
        curvature=summary['curvature'],
        rv=summary['rv'],
        mass_pi1=summary['mass_pi1'],
        rms=summary['rms'],
        rms_all=summary['rms_all'],
        n_obs=summary['n_obs'],
        n_excluded=summary['n_excluded'],
        residuals=residuals,
    )


def bootstrap_orbit_elements(result, n_boot=100, rng=None, verbose=False, progress=True):
    """
    Parametric-bootstrap uncertainty estimate for a fitted orbit.

    Takes the fitted model's predicted sky position (rho_c, theta_c) at each
    observation epoch from result.residuals, converts it to tangent-plane
    (x, y) coordinates, and draws n_boot synthetic observation sets by adding
    isotropic Gaussian noise with std = result.rms (arcsec) to (x, y) —
    matching the noise model the RMS itself is defined against. Each
    synthetic set is converted back to (rho, theta) and refit from scratch
    with the same procedure as fit_orbit_dataframe (full global search +
    outlier rejection), at the original observation epochs.

    result   — an OrbitFitResult from fit_orbit_dataframe(...).
    n_boot   — number of bootstrap realizations (default 100).
    rng      — optional numpy.random.Generator, for reproducibility.
    verbose  — print one line per realization (routed through the progress
               bar's own write, if one is active, so it doesn't clobber it).
    progress — show a live progress bar (tqdm.auto — a plain bar in a
               terminal, a widget in Jupyter). Silently falls back to no bar
               if the optional 'tqdm' package isn't installed.

    Returns a DataFrame with one row per successful bootstrap fit and
    columns 'a', 'e', 'T', 'P', 'omega', 'Node', 'i', 'rms'. Realizations
    where the refit raises (e.g. the global search finds no valid orbit)
    are silently skipped — n_boot is an upper bound on the row count, not
    a guarantee.

    Note: this reruns the full fit_orbit_dataframe pipeline n_boot times,
    so it costs roughly n_boot times a single fit (seconds to minutes for
    n_boot=100, depending on the number of observations).
    """
    if rng is None:
        rng = np.random.default_rng()

    bar = None
    if progress:
        try:
            from tqdm.auto import tqdm
            bar = tqdm(total=n_boot, desc='bootstrap orbits')
        except ImportError:
            pass

    def log(msg):
        bar.write(msg) if bar is not None else print(msg)

    epochs_abs = result.residuals['epoch'].to_numpy(dtype=float)
    rho_c      = result.residuals['rho_c'].to_numpy(dtype=float)
    theta_c    = result.residuals['theta_c'].to_numpy(dtype=float)
    x_c = rho_c * np.sin(np.radians(theta_c))
    y_c = rho_c * np.cos(np.radians(theta_c))
    sigma = result.rms

    rows = []
    try:
        for b in range(n_boot):
            x_sim = x_c + rng.normal(0.0, sigma, size=x_c.shape)
            y_sim = y_c + rng.normal(0.0, sigma, size=y_c.shape)
            rho_sim   = np.hypot(x_sim, y_sim)
            theta_sim = np.degrees(np.arctan2(x_sim, y_sim)) % 360.0

            boot_df = pd.DataFrame({'year': epochs_abs, 'rho': rho_sim, 'theta': theta_sim})
            try:
                boot_result = fit_orbit_dataframe(boot_df, year_col='year', rho_col='rho',
                                                  theta_col='theta', verbose=False)
            except Exception as exc:
                if verbose:
                    log(f"  bootstrap {b+1}/{n_boot}: FAILED ({exc})")
                if bar is not None:
                    bar.update(1)
                continue

            row = dict(boot_result.elements)
            row['rms'] = boot_result.rms
            rows.append(row)
            if verbose:
                log(f"  bootstrap {b+1}/{n_boot}: a={row['a']:.4f}  e={row['e']:.4f}  "
                   f"P={row['P']:.3f}  rms={row['rms']:.5f}")
            if bar is not None:
                bar.set_postfix(ok=len(rows), a=f"{row['a']:.3g}", P=f"{row['P']:.3g}")
                bar.update(1)
    finally:
        if bar is not None:
            bar.close()

    return pd.DataFrame(rows, columns=list(ORBIT_ELEMENTS) + ['rms'])


def orbit_element_limits(bootstrap_df, elements=ORBIT_ELEMENTS, lower_pct=16.0, upper_pct=84.0):
    """
    Percentile-based upper/lower limits for each orbital element from a
    bootstrap sample (see bootstrap_orbit_elements).

    lower_pct, upper_pct — default 16/84, i.e. the ~1-sigma interval for a
    roughly Gaussian bootstrap distribution. Pass (0, 100) for the full
    range, or (2.5, 97.5) for a ~95% interval.

    Returns a DataFrame indexed by element name with columns 'lower' and
    'upper'.

    Note: omega and Node are angles (deg, 0-360) and can wrap around 0/360.
    If a bootstrap distribution straddles the wrap point, these plain
    percentiles will be misleading — check the corner plot in that case.
    """
    limits = {}
    for el in elements:
        vals = bootstrap_df[el].dropna()
        limits[el] = {
            'lower': np.percentile(vals, lower_pct),
            'upper': np.percentile(vals, upper_pct),
        }
    return pd.DataFrame(limits).T[['lower', 'upper']]


def plot_orbit_corner(bootstrap_df, elements=ORBIT_ELEMENTS, best_fit=None,
                      filename='orbit_corner.png', dpi=150, log_elements=('a', 'P')):
    """
    Save a corner plot (pairwise scatter + 1-D histograms) of the bootstrap
    orbital-element distribution to a PNG file.

    bootstrap_df  — DataFrame from bootstrap_orbit_elements(...).
    elements      — which columns to plot, and their order.
    best_fit      — optional dict of the nominal (non-bootstrap) element
                    values, e.g. result.elements — drawn as red reference
                    lines on every panel.
    filename      — output PNG path.
    log_elements  — element names to plot on a log scale (default: a, P —
                    both strictly positive and, for short observation arcs,
                    often degenerate across orders of magnitude, which
                    crushes a linear-axis panel into one corner).

    Returns filename.
    """
    import matplotlib.pyplot as plt

    SURFACE  = '#fcfcfb'
    INK      = '#0b0b0b'
    INK_SOFT = '#52514e'
    MUTED    = '#898781'
    GRID     = '#e1e0d9'
    SAMPLE   = '#2a78d6'   # bootstrap sample (categorical slot 1, blue)
    NOMINAL  = '#e34948'   # best (non-bootstrap) fit (categorical slot 8, red)

    n = len(elements)
    fig, axes = plt.subplots(n, n, figsize=(2.3 * n, 2.3 * n), facecolor=SURFACE)

    for i in range(n):
        yi = elements[i]
        for j in range(n):
            xi = elements[j]
            ax = axes[i, j]
            ax.set_facecolor(SURFACE)
            for spine in ax.spines.values():
                spine.set_color(GRID)

            if j > i:
                ax.axis('off')
                continue

            x_log = xi in log_elements
            y_log = yi in log_elements

            if i == j:
                vals = bootstrap_df[xi].dropna()
                if x_log:
                    bins = np.geomspace(vals.min(), vals.max(), 21)
                    ax.set_xscale('log')
                else:
                    bins = 20
                ax.hist(vals, bins=bins, color=SAMPLE, alpha=0.75, edgecolor=SURFACE)
                lo, hi = np.percentile(vals, [16, 84])
                for v in (lo, hi):
                    ax.axvline(v, color=MUTED, lw=1, ls='--')
                if best_fit is not None and xi in best_fit:
                    ax.axvline(best_fit[xi], color=NOMINAL, lw=1.5)
                ax.set_yticks([])
            else:
                if x_log:
                    ax.set_xscale('log')
                if y_log:
                    ax.set_yscale('log')
                ax.scatter(bootstrap_df[xi], bootstrap_df[yi], s=8, color=SAMPLE,
                          alpha=0.35, edgecolors='none')
                if best_fit is not None and xi in best_fit and yi in best_fit:
                    ax.axvline(best_fit[xi], color=NOMINAL, lw=1, alpha=0.6)
                    ax.axhline(best_fit[yi], color=NOMINAL, lw=1, alpha=0.6)
                ax.grid(True, color=GRID, lw=0.5)

            ax.tick_params(colors=MUTED, labelsize=7)
            if i < n - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(ELEMENT_LABELS.get(xi, xi), color=INK_SOFT, fontsize=9)
            if j > 0 or i == j:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(ELEMENT_LABELS.get(yi, yi), color=INK_SOFT, fontsize=9)

    fig.suptitle(f'Orbital-element bootstrap distribution (N={len(bootstrap_df)})',
                color=INK, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(filename, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)
    return filename


def plot_orbit_family(bootstrap_df, result, filename='orbit_family.png',
                      n_orbits=None, n_points=2000, dpi=150, zoom_pad=0.15,
                      point_size=1.2):
    """
    Save a plot of the fitted orbit family in the plane of the sky (tangent
    plane, arcsec) to a PNG file: each bootstrap realization in bootstrap_df
    drawn as a scatter of points sampled around its orbit (not a connected
    curve — for a highly eccentric/degenerate orbit, evenly time-sampled
    points are far apart near periastron, and connecting them with straight
    lines draws a jagged artifact, especially where the curve is only
    partly inside the view). The nominal fit is still drawn as a curve
    (it's a single well-behaved orbit), and the observations (from
    result.residuals) as points — hollow for any flagged excluded=True.

    bootstrap_df — DataFrame from bootstrap_orbit_elements(...); each row's
                   Campbell elements are turned back into an orbit curve via
                   tilin.campbell_to_ti.
    result       — the OrbitFitResult the bootstrap was built from (supplies
                   the nominal elements and the observed positions).
    n_orbits     — plot at most this many bootstrap orbits (default: all).
    n_points     — points sampled per orbit, over one full period. Points
                   outside the view window are simply not visible — no
                   orbit is dropped or specially handled, unlike a
                   connected-line rendering.
    point_size   — matplotlib scatter marker size for the sampled points.
    filename     — output PNG path.
    zoom_pad     — fractional padding for the view window, sized to the
                   nominal orbit curve + observations. Pass None to instead
                   autoscale to the full bootstrap family.

    Returns filename.
    """
    import matplotlib.pyplot as plt

    SURFACE  = '#fcfcfb'
    INK      = '#0b0b0b'
    INK_SOFT = '#52514e'
    MUTED    = '#898781'
    GRID     = '#e1e0d9'
    SAMPLE   = '#2a78d6'   # bootstrap orbit family (categorical slot 1, blue)
    NOMINAL  = '#e34948'   # best (non-bootstrap) fit (categorical slot 8, red)

    sample_rows = bootstrap_df
    if n_orbits is not None:
        sample_rows = sample_rows.head(n_orbits)

    fig, ax = plt.subplots(figsize=(7, 7), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_color(GRID)

    def orbit_xy(a, e, T, P, omega, Node, i):
        ti = _core.campbell_to_ti(a, e, T, P, omega, Node, i)
        t  = np.linspace(T - P / 2.0, T + P / 2.0, n_points)
        rho   = np.empty(n_points)
        theta = np.empty(n_points)
        for k, tk in enumerate(t):
            rho[k], theta[k] = _core.rtti(ti, tk)
        x = rho * np.sin(np.radians(theta))
        y = rho * np.cos(np.radians(theta))
        return x, y

    for _, row in sample_rows.iterrows():
        x, y = orbit_xy(row['a'], row['e'], row['T'], row['P'], row['omega'], row['Node'], row['i'])
        ax.scatter(x, y, s=point_size, color=SAMPLE, alpha=0.25, edgecolors='none')

    el = result.elements
    x_nom, y_nom = orbit_xy(el['a'], el['e'], el['T'], el['P'], el['omega'], el['Node'], el['i'])
    ax.plot(x_nom, y_nom, color=NOMINAL, lw=1.8, label='nominal fit')

    used = result.residuals[~result.residuals['excluded']]
    excl = result.residuals[result.residuals['excluded']]
    ax.scatter(used['x'], used['y'], s=28, color=INK, zorder=5, label='observed')
    if len(excl):
        ax.scatter(excl['x'], excl['y'], s=28, facecolors='none', edgecolors=MUTED,
                  zorder=5, label='observed (excluded)')

    ax.plot(0, 0, marker='+', color=INK, ms=10, mew=1.5, zorder=6)

    if zoom_pad is not None:
        view_x = np.concatenate([x_nom, used['x'].to_numpy(), excl['x'].to_numpy()])
        view_y = np.concatenate([y_nom, used['y'].to_numpy(), excl['y'].to_numpy()])
        pad_x = zoom_pad * (view_x.max() - view_x.min()) or 1.0
        pad_y = zoom_pad * (view_y.max() - view_y.min()) or 1.0
        ax.set_xlim(view_x.min() - pad_x, view_x.max() + pad_x)
        ax.set_ylim(view_y.min() - pad_y, view_y.max() + pad_y)

    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color=GRID, lw=0.5)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.set_xlabel('X (arcsec)', color=INK_SOFT, fontsize=10)
    ax.set_ylabel('Y (arcsec)', color=INK_SOFT, fontsize=10)
    ax.set_title(f'Orbit family in the plane of the sky (N={len(sample_rows)} bootstrap orbits)',
                color=INK, fontsize=11)
    legend = ax.legend(fontsize=8, facecolor=SURFACE, edgecolor=GRID, labelcolor=INK_SOFT)

    fig.tight_layout()
    fig.savefig(filename, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)
    return filename
