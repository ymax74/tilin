#!/usr/bin/env python3
"""
TILIN132 — Thiele-Innes method for visual binary star orbit fitting.
Translated from the Pulkovo Observatory Pascal program TILIN132.

Input:  input.txt  — columns: epoch(yr)  theta(deg)  rho(arcsec)
Output: output.txt — fitted Campbell elements + residuals table

Usage:
    python tilin.py
(run from the directory that contains input.txt)
"""

import sys
import math
import numpy as np
from scipy.linalg import lstsq
from scipy.optimize import least_squares

# ---------------------------------------------------------------------------
# Conversion constants
# ---------------------------------------------------------------------------
DEG = math.pi / 180.0
# AU/yr speed factor for radial velocity conversion (arcsec/yr → AU/yr)
YEAR_TO_AU = 3600.0 * 24.0 * 365.2422 / 149600000.0

# ---------------------------------------------------------------------------
# Angle utilities
# ---------------------------------------------------------------------------

def atan2_pa(x, y):
    """Position-angle atan2: arctan(x/y) with quadrant from sign of y, result in [0, 2π)."""
    a = math.atan2(x, y)
    if a < 0.0:
        a += 2.0 * math.pi
    return a


def norm_deg(a):
    """Reduce angle to (-180, 180] degrees — used for residual wrapping."""
    while a > 180.0:
        a -= 360.0
    while a <= -180.0:
        a += 360.0
    return a


def norm_rad(a):
    """Reduce angle to [0, 2π)."""
    pi2 = 2.0 * math.pi
    while a >= pi2:
        a -= pi2
    while a < 0.0:
        a += pi2
    return a


# ---------------------------------------------------------------------------
# Kepler equation
# ---------------------------------------------------------------------------

def solve_kepler(M, e, eps=1e-13):
    """Solve E - e·sin(E) = M  (elliptic, e<1) or  e·sinh(H) - H = M (hyperbolic)."""
    if e >= 1.0:
        H = M
        for _ in range(50):
            fh = e * math.sinh(H) - H - M
            dH = fh / (e * math.cosh(H) - 1.0)
            H -= dH
            if abs(dH) < eps:
                break
        return H
    else:
        pi = math.pi
        while M > pi:
            M -= 2.0 * pi
        while M < -pi:
            M += 2.0 * pi
        E = M
        for _ in range(50):
            fe = E - e * math.sin(E) - M
            dE = fe / (1.0 - e * math.cos(E))
            E -= dE
            if abs(dE) < eps:
                break
        return E


# ---------------------------------------------------------------------------
# Predict (rho, theta) from orbital state vector ti
#
#  ti[0] = abp   — semi-major axis (arcsec)
#  ti[1] = e     — eccentricity
#  ti[2] = th    — time of periastron relative to t0 (years)
#  ti[3] = prd   — period (years)
#  ti[4] = aa  \
#  ti[5] = bb   | Thiele-Innes direction cosines (unit-free; abp absorbed)
#  ti[6] = ff   |
#  ti[7] = gg  /
#  ti[8] = omm  — argument of periastron (rad)
#  ti[9] = omb  — longitude of ascending node (rad)
#  ti[10]= ai   — inclination (rad)
# ---------------------------------------------------------------------------

def rtti(ti, t):
    """Given state vector ti and epoch t (relative to t0), return (rho, theta_deg)."""
    abp, e, th, prd, aa, bb, ff, gg = ti[0], ti[1], ti[2], ti[3], ti[4], ti[5], ti[6], ti[7]
    M = (t - th) * 2.0 * math.pi / prd
    E = solve_kepler(M, e)
    if e < 1.0:
        sE = math.sin(E);  cE = math.cos(E)
        xi = cE - e
    else:
        sE = math.sinh(E); cE = math.cosh(E)
        xi = -cE + e
    yi = math.sqrt(abs(1.0 - e * e)) * sE
    xc = abs(abp) * (aa * xi + ff * yi)
    yc = abs(abp) * (bb * xi + gg * yi)
    rho   = math.sqrt(xc * xc + yc * yc)
    theta = atan2_pa(xc, yc) * 180.0 / math.pi
    return rho, theta


# ---------------------------------------------------------------------------
# 3-D Cartesian state → Keplerian elements  (elor3 / elor2 logic)
# ---------------------------------------------------------------------------

def cartesian_to_kepler(x, y, z, vx, vy, vz, kp):
    """
    Convert Cartesian state (AU, AU/yr) to Keplerian elements.
    kp = sqrt(4π²·mass) in AU^(3/2)/yr units.
    Returns (a_AU, e, omm, omb, ai, th) or None on failure.
    omb is adjusted as in Pascal: omb = π/2 - omb_raw.
    """
    r  = math.sqrt(x*x + y*y + z*z)
    v2 = vx*vx + vy*vy + vz*vz
    rdv = x*vx + y*vy + z*vz

    # Semi-major axis
    inv_a = 2.0/r - v2/(kp*kp)
    if abs(inv_a) < 1e-30:
        return None
    a = 1.0 / inv_a

    pf = v2*r/(kp*kp) - 1.0       # (a - r)/a  = 1 - r/a for elliptic
    denom = kp*kp * abs(a)
    if denom < 1e-30:
        return None
    e2 = rdv*rdv / denom + pf*pf
    e  = math.sqrt(abs(e2))

    # Eccentric anomaly at epoch
    if e > 1e-10:
        sE0 = rdv / (kp * math.sqrt(abs(a)) * e)
        cE0 = pf / e
    else:
        sE0 = 0.0;  cE0 = 1.0

    # Unit vectors along periastron (px,py,pz) and 90° ahead (qx,qy,qz)
    ak = cE0 / r
    al = math.sqrt(abs(a)) * sE0 / kp
    px = ak*x - al*vx
    py = ak*y - al*vy
    pz = ak*z - al*vz

    if e > 1e-10:
        fac_q = r * math.sqrt(abs(1.0 - e*e))
        if abs(fac_q) < 1e-30:
            return None
        am = sE0 / fac_q
        if a > 0.0:
            an = math.sqrt(abs(a / (1.0 - e*e))) * (cE0 - e) / kp
        else:
            an = -math.sqrt(abs(a / (1.0 - e*e))) * (cE0 - e) / kp
    else:
        am = 0.0;  an = 0.0

    qx = am*x + an*vx
    qy = am*y + an*vy
    qz = am*z + an*vz

    # Inclination from pz, qz
    si = math.sqrt(pz*pz + qz*qz)
    if si > 1e-10:
        somm = pz / si
        comm = qz / si
        omm  = math.asin(max(-1.0, min(1.0, somm)))
        if comm < 0.0:
            omm = math.pi - omm
        if omm < 0.0:
            omm += 2.0*math.pi

        somb = py*comm - qy*somm
        comb = px*comm - qx*somm
        omb  = math.asin(max(-1.0, min(1.0, somb)))
        if comb < 0.0:
            omb = math.pi - omb
        if omb < 0.0:
            omb += 2.0*math.pi

        ai = math.asin(max(-1.0, min(1.0, si)))
        if x*vy - y*vx < 0.0:
            ai = math.pi - ai
    else:
        omm = 0.5;  omb = 0.5;  ai = 0.5

    # Epoch of periastron (relative)
    if e < 1.0:
        E0 = math.asin(max(-1.0, min(1.0, sE0)))
        if cE0 < 0.0:
            E0 = math.pi - E0
        if E0 < 0.0:
            E0 += 2.0*math.pi
        n  = kp / (abs(a) ** 1.5)          # mean motion
        th = -(E0 - e*sE0) / n
    else:
        H0 = math.log(sE0 + math.sqrt(sE0*sE0 + 1.0)) if abs(sE0) < 1e15 else sE0
        n  = kp / (abs(a) ** 1.5)
        th = (H0 - e*sE0) / n              # hyperbolic: sign flipped

    # Apply Pascal convention: omb → π/2 - omb
    omb = norm_rad(math.pi/2.0 - omb)

    return abs(a), e, omm, omb, ai, th


# ---------------------------------------------------------------------------
# Phase-space state → orbital state vector ti  (tpvd equivalent)
# ---------------------------------------------------------------------------

def tpvd(rho0, theta0, rho_dot, theta_dot, curvature, vr_arcsec, mass, t_ref, parallax):
    """
    Convert initial-condition phase-space state to orbital state vector ti.

    rho0       — separation at t_ref (arcsec)
    theta0     — position angle at t_ref (degrees)
    rho_dot    — dρ/dt (arcsec/yr)
    theta_dot  — dθ/dt (deg/yr)
    curvature  — radius of curvature of apparent path (arcsec)
    vr_arcsec  — radial velocity (arcsec/yr)
    mass       — total mass (solar units, for parallax=1 AU)
    t_ref      — not used (kept for API symmetry)
    parallax   — parallax (arcsec); 1 means distance = 1 pc

    Returns ti[0..10] or None on failure.
    """
    t0r = theta0 * DEG

    # Sky-plane Cartesian position (arcsec)
    xs_sky = rho0 * math.sin(t0r)
    ys_sky = rho0 * math.cos(t0r)

    # Sky-plane velocity (arcsec/yr)
    vx_sky = rho_dot * math.sin(t0r) + rho0 * theta_dot * math.cos(t0r) * DEG
    vy_sky = rho_dot * math.cos(t0r) - rho0 * theta_dot * math.sin(t0r) * DEG

    mu  = math.sqrt(vx_sky**2 + vy_sky**2)     # speed in sky plane (arcsec/yr)
    psi = atan2_pa(vx_sky, vy_sky)              # direction of velocity (rad)

    kp  = math.sqrt(4.0 * math.pi**2 * mass)

    # Convert to AU
    xa = xs_sky / parallax
    ya = ys_sky / parallax

    # 3-D distance from radius of curvature
    sin_ang = abs(math.sin(psi - t0r))
    if mu > 1e-30 and curvature > 1e-10 and sin_ang > 1e-30:
        arg = kp**2 * rho0 * curvature * sin_ang / mu**2
        rt  = abs(arg) ** (1.0/3.0)
    else:
        rt = rho0 / parallax    # fallback: put it at the 2-D distance

    ro_au  = rho0 / parallax
    za_sq  = rt**2 - ro_au**2
    za     = math.sqrt(abs(za_sq))              # z-component (AU), sign = +

    # 3-D velocity (AU/yr)
    xt_au = mu * math.sin(psi) / parallax
    yt_au = mu * math.cos(psi) / parallax
    zt_au = vr_arcsec * YEAR_TO_AU             # radial velocity

    result = cartesian_to_kepler(xa, ya, za, xt_au, yt_au, zt_au, kp)
    if result is None:
        return None

    a_au, e, omm, omb, ai, th = result

    if a_au <= 1e-15 or e > 10.0:
        return None

    # Period (years)
    prd = 2.0 * math.pi * a_au**1.5 / kp

    if prd < 1.0:
        return None

    # Semi-major axis in arcsec
    abp = a_au * parallax

    # Thiele-Innes direction cosine matrix
    # NOTE: Pascal convention: comb = sin(omb), somb = cos(omb)
    comm = math.cos(omm);  somm = math.sin(omm)
    comb = math.sin(omb);  somb = math.cos(omb)   # Pascal naming
    ci   = math.cos(ai)

    aa = comm*comb - somm*somb*ci
    bb = comm*somb + somm*comb*ci
    ff = -somm*comb - comm*somb*ci
    gg = -somm*somb + comm*comb*ci

    ti = [abp, e, th, prd, aa, bb, ff, gg, omm, omb, ai]
    return ti


# ---------------------------------------------------------------------------
# Campbell elements from ABFG direction cosines (GetElFrTi)
# ---------------------------------------------------------------------------

def angles_from_abfg_cosines(aa, bb, ff, gg):
    """
    Recover (omm, omb, ai) in radians from unit-vector Thiele-Innes coefficients.
    This implements GetElFrTi from MPDQ.PAS.
    """
    rr = 0.5 * (aa*aa + bb*bb + ff*ff + gg*gg)
    ss = aa*gg - bb*ff
    inner = (rr + ss) * (rr - ss)
    abp_sq = rr + math.sqrt(max(0.0, inner))     # = abp_unit²  (abp_unit = 1 here)

    if abp_sq < 1e-20:
        return 0.0, 0.0, 0.0

    cosai = (aa*gg - bb*ff) / abp_sq
    cosai = max(-1.0, min(1.0, cosai))
    ai = math.acos(cosai)

    ombMomm = atan2_pa(bb + ff, aa - gg)
    ombPomm = atan2_pa(bb - ff, aa + gg)

    omb = (ombPomm + ombMomm) * 0.5
    omb = norm_rad(omb)
    omm = norm_rad(ombPomm - omb)

    # Pascal convention: omb → π/2 - omb
    omb = norm_rad(math.pi/2.0 - omb)

    return omm, omb, ai


def campbell_from_ti(ti, t0):
    """Extract Campbell elements dict from state vector ti and reference epoch t0."""
    abp = ti[0];  e = ti[1];  th = ti[2];  prd = ti[3]
    omm = ti[8];  omb = ti[9];  ai = ti[10]

    # Inclination as reported in Pascal: 180 - ai_deg
    i_deg = 180.0 - ai * 180.0 / math.pi
    if i_deg < 0.0:
        i_deg += 360.0

    return {
        'a':     abp,
        'e':     e,
        'T':     th + t0,
        'P':     prd,
        'omega': omm * 180.0 / math.pi,
        'Node':  omb * 180.0 / math.pi,
        'i':     i_deg,
    }


def campbell_to_ti(a, e, T, P, omega_deg, Node_deg, i_deg):
    """
    Build a state vector ti from Campbell elements — the inverse of
    campbell_from_ti (given the same t0/T reference). Only (t - T) enters
    rtti(ti, t), so T and the epochs passed to rtti may be either both
    absolute or both relative to the same t0.
    """
    omm = omega_deg * DEG
    omb = Node_deg * DEG
    ai  = (180.0 - i_deg) * DEG   # inverse of the i_deg = 180 - ai*180/pi convention

    comm = math.cos(omm);  somm = math.sin(omm)
    comb = math.sin(omb);  somb = math.cos(omb)   # Pascal naming convention, as in tpvd
    ci   = math.cos(ai)

    aa =  comm*comb - somm*somb*ci
    bb =  comm*somb + somm*comb*ci
    ff = -somm*comb - comm*somb*ci
    gg = -somm*somb + comm*comb*ci

    return [a, e, T, P, aa, bb, ff, gg, omm, omb, ai]


# ---------------------------------------------------------------------------
# Least-squares fit for Thiele-Innes scaled constants (A, B, F, G = abp·aa, ...)
# given fixed (e, th, prd)   — GetPrmElOsn from MPDQ.PAS
# ---------------------------------------------------------------------------

def fit_thiele_innes(epochs, rho_obs, theta_obs, weights, excluded, e, th, prd):
    """
    For fixed (e, th, prd) solve the linear system for A=abp·aa, B=abp·bb, F=abp·ff, G=abp·gg.
    Returns (A, B, F, G, rms_arcsec) or None.

    Model: X_i = A·xi + F·yi,  Y_i = B·xi + G·yi
    where xi = cos(E)-e,  yi = sqrt(|1-e²|)·sin(E)  from Kepler at epoch t_i.
    """
    pi2prd = 2.0 * math.pi / prd
    Mx_rows = [];  bx = [];  My_rows = [];  by_arr = []
    for k in range(len(epochs)):
        if excluded[k]:
            continue
        t  = epochs[k];  rho = rho_obs[k];  theta = theta_obs[k]
        w  = weights[k]
        M  = (t - th) * pi2prd
        E  = solve_kepler(M, e)
        if e < 1.0:
            xi = math.cos(E) - e
            yi = math.sqrt(abs(1.0 - e*e)) * math.sin(E)
        else:
            xi = -(math.cosh(E) - e)
            yi = math.sqrt(abs(e*e - 1.0)) * math.sinh(E)

        X_obs = rho * math.sin(theta * DEG)
        Y_obs = rho * math.cos(theta * DEG)

        Mx_rows.append([xi*w, yi*w])
        bx.append(X_obs * w)
        My_rows.append([xi*w, yi*w])
        by_arr.append(Y_obs * w)

    n = len(bx)
    if n < 2:
        return None

    Mx  = np.array(Mx_rows);  bxv = np.array(bx)
    My  = np.array(My_rows);  byv = np.array(by_arr)

    sol_x, _, _, _ = lstsq(Mx, bxv)   # [A, F]
    sol_y, _, _, _ = lstsq(My, byv)   # [B, G]

    A, F = float(sol_x[0]), float(sol_x[1])
    B, G = float(sol_y[0]), float(sol_y[1])

    # RMS residual in XY plane
    res2 = 0.0
    nn   = 0
    pi2prd = 2.0 * math.pi / prd
    for k in range(len(epochs)):
        if excluded[k]:
            continue
        t   = epochs[k];  rho = rho_obs[k];  theta = theta_obs[k]
        M   = (t - th) * pi2prd
        E   = solve_kepler(M, e)
        if e < 1.0:
            xi = math.cos(E) - e
            yi = math.sqrt(abs(1.0 - e*e)) * math.sin(E)
        else:
            xi = -(math.cosh(E) - e)
            yi = math.sqrt(abs(e*e - 1.0)) * math.sinh(E)
        xc = A*xi + F*yi;  yc = B*xi + G*yi
        X_obs = rho * math.sin(theta * DEG)
        Y_obs = rho * math.cos(theta * DEG)
        res2 += (X_obs - xc)**2 + (Y_obs - yc)**2
        nn += 1

    oev = math.sqrt(res2 / (2*nn)) if nn > 0 else 1e38
    return A, B, F, G, oev


# ---------------------------------------------------------------------------
# Build state vector ti from scaled Thiele-Innes constants A=abp·aa, …
# ---------------------------------------------------------------------------

def ti_from_scaled(e, th, prd, A, B, F, G):
    """Build ti state vector from scaled Thiele-Innes constants."""
    rr  = 0.5 * (A*A + B*B + F*F + G*G)
    ss  = A*G - B*F
    inner = (rr + ss) * (rr - ss)
    abp = math.sqrt(rr + math.sqrt(max(0.0, inner)))
    if abp < 1e-15:
        return None
    aa = A/abp;  bb = B/abp;  ff = F/abp;  gg = G/abp
    omm, omb, ai = angles_from_abfg_cosines(aa, bb, ff, gg)
    return [abp, e, th, prd, aa, bb, ff, gg, omm, omb, ai]


# ---------------------------------------------------------------------------
# Weighted RMS residual (fnev equivalent)
# ---------------------------------------------------------------------------

def rms_residual(ti, epochs, rho_obs, theta_obs, weights, excluded):
    """Weighted RMS matching Pascal fnev: sqrt(0.5·Σw²(δρ²+δθ_arc²)/n_used)."""
    total = 0.0;  n = 0
    for k in range(len(epochs)):
        if excluded[k]:
            continue
        try:
            rho_c, theta_c = rtti(ti, epochs[k])
        except Exception:
            return 1e38
        drho   = rho_c - rho_obs[k]
        dtheta = norm_deg(theta_c - theta_obs[k])
        dte    = dtheta * math.pi * rho_obs[k] / 180.0
        w      = weights[k]
        total += w*w * (drho*drho + dte*dte)
        n     += 1
    if n == 0:
        return 1e38
    return math.sqrt(0.5 * total / n)


def rms_all_obs(ti, epochs, rho_obs, theta_obs, weights):
    """RMS over all observations (for 'all obs.' output line)."""
    total = 0.0;  n = len(epochs)
    for k in range(n):
        try:
            rho_c, theta_c = rtti(ti, epochs[k])
        except Exception:
            return 1e38
        drho   = rho_obs[k] - rho_c
        dtheta = norm_deg(theta_obs[k] - theta_c)
        dte    = dtheta * math.pi * rho_obs[k] / 180.0
        w      = weights[k]
        total += w*w * (drho*drho + dte*dte)
    return math.sqrt(0.5 * total / n) if n > 0 else 1e38


# ---------------------------------------------------------------------------
# Residual vector for scipy.optimize.least_squares
# params = [e, th, prd, A, B, F, G]
# ---------------------------------------------------------------------------

def residuals_for_lsq(params, epochs, rho_obs, theta_obs, weights, excluded):
    e, th, prd, A, B, F, G = params
    if prd < 0.5 or e < 0.0 or e >= 1.0:
        n = sum(1 for ex in excluded if not ex)
        return np.ones(2 * max(n, 1)) * 1e6

    ti = ti_from_scaled(e, th, prd, A, B, F, G)
    if ti is None:
        n = sum(1 for ex in excluded if not ex)
        return np.ones(2 * max(n, 1)) * 1e6

    resids = []
    for k in range(len(epochs)):
        if excluded[k]:
            continue
        try:
            rho_c, theta_c = rtti(ti, epochs[k])
        except Exception:
            resids.extend([1e6, 1e6])
            continue
        drho   = rho_c - rho_obs[k]
        dtheta = norm_deg(theta_c - theta_obs[k])
        dte    = dtheta * math.pi * rho_obs[k] / 180.0
        w      = weights[k]
        resids.extend([w * drho, w * dte])

    return np.array(resids) if resids else np.array([0.0])


# ---------------------------------------------------------------------------
# Local refinement with scipy.optimize.least_squares (Levenberg-Marquardt)
# ---------------------------------------------------------------------------

def local_refine(init_params, epochs, rho_obs, theta_obs, weights, excluded,
                 max_nfev=50000, tol=1e-11):
    """
    Refine (e, th, prd, A, B, F, G) using trust-region reflective LSQ.
    Returns (params, ti, rms) or (None, None, 1e38) on failure.
    """
    x0 = np.array(init_params, dtype=float)
    bounds_lo = [0.0,   -np.inf, 0.5,   -np.inf, -np.inf, -np.inf, -np.inf]
    bounds_hi = [0.9999, np.inf, np.inf,  np.inf,  np.inf,  np.inf,  np.inf]

    try:
        result = least_squares(
            residuals_for_lsq, x0,
            args=(epochs, rho_obs, theta_obs, weights, excluded),
            bounds=(bounds_lo, bounds_hi),
            method='trf',
            max_nfev=max_nfev,
            ftol=tol, xtol=tol, gtol=tol,
            verbose=0,
        )
        params = result.x.tolist()
    except Exception:
        return None, None, 1e38

    ti  = ti_from_scaled(*params)
    if ti is None:
        return None, None, 1e38
    rms = rms_residual(ti, epochs, rho_obs, theta_obs, weights, excluded)
    return params, ti, rms


# ---------------------------------------------------------------------------
# Outlier rejection (maxiskl4-style IQR-based sigma)
# ---------------------------------------------------------------------------

OUTLIER_SIGMA = 5.0    # Pascal: porog = 5

def find_worst_outlier(ti, epochs, rho_obs, theta_obs, excluded):
    """
    Find the observation with the largest normalised residual (IQR method).
    Returns (index, max_sigma) or (-1, 0) if nothing to flag.
    """
    drho_list  = [];  dtheta_list = []
    indices    = []
    for k in range(len(epochs)):
        if excluded[k]:
            continue
        try:
            rho_c, theta_c = rtti(ti, epochs[k])
        except Exception:
            continue
        drho_list.append(rho_c - rho_obs[k])
        dtheta_list.append(norm_deg(theta_c - theta_obs[k]))
        indices.append(k)

    if len(drho_list) < 4:
        return -1, 0.0

    # IQR-based σ using the 1-σ percentiles (intver = (1-2*0.34134)/2 ≈ 0.1586)
    intver = (1.0 - 2.0 * 0.34134) / 2.0
    n = len(drho_list)
    p1 = max(0, min(n-1, int(n * intver)))
    p2 = max(0, min(n-1, int(n * (1.0 - intver))))

    drho_sorted   = sorted(drho_list)
    dtheta_sorted = sorted(dtheta_list)

    sigma_rho   = (drho_sorted[p2]   - drho_sorted[p1])   * 0.5
    sigma_theta = (dtheta_sorted[p2] - dtheta_sorted[p1]) * 0.5

    if sigma_rho   < 1e-10: sigma_rho   = 1e-10
    if sigma_theta < 1e-10: sigma_theta = 1e-10

    best_k = -1;  best_nev = 0.0
    for j, k in enumerate(indices):
        nev_r = abs(drho_list[j]) / sigma_rho
        nev_t = abs(dtheta_list[j]) / sigma_theta
        nev   = max(nev_r, nev_t)
        if nev > best_nev:
            best_nev = nev;  best_k = k

    return best_k, best_nev


# ---------------------------------------------------------------------------
# Pascal-faithful global search  (Pascal: odflM4d in DERTY2.PAS)
#
# The real TILIN132 program does NOT use the circular-arc/ciriteracii path
# below for a generic input (mbcir=false) — it draws 2000 random trial
# orbits, scores each with GetPrmEl (== fit_thiele_innes), refines the best
# 10 with a Newton/LM differential correction (elemiter2, substituted here
# by local_refine — both converge to the same stationary point given the
# same starting candidate), and keeps the lowest-RMS result.  This runs
# twice in a row (two odflM4d calls in TILIN132.PAS's main) sharing one
# continuous RNG stream; the *second* pass's winner is what gets written
# out.  random() itself is Free Pascal's MT19937 (verified empirically
# against a locally-compiled TILIN132 binary: RandSeed=0, classic
# init_genrand seeding, one raw tempered 32-bit word per call / 2**32).
# ---------------------------------------------------------------------------

class PascalRandom:
    """Free Pascal's `random` (extended/real form): MT19937, RandSeed=0,
    classic init_genrand seeding, one tempered 32-bit word per call."""

    def __init__(self, seed=0):
        self.mt = [0] * 624
        self.index = 624
        self.mt[0] = seed & 0xFFFFFFFF
        for i in range(1, 624):
            self.mt[i] = (1812433253 * (self.mt[i-1] ^ (self.mt[i-1] >> 30)) + i) & 0xFFFFFFFF

    def _generate(self):
        for i in range(624):
            y = (self.mt[i] & 0x80000000) + (self.mt[(i+1) % 624] & 0x7fffffff)
            self.mt[i] = self.mt[(i + 397) % 624] ^ (y >> 1)
            if y % 2 != 0:
                self.mt[i] ^= 2567483615
        self.index = 0

    def _next_word(self):
        if self.index >= 624:
            self._generate()
        y = self.mt[self.index]
        y ^= y >> 11
        y ^= (y << 7) & 2636928640
        y ^= (y << 15) & 4022730752
        y ^= y >> 18
        self.index += 1
        return y

    def random(self):
        """Pascal `random` (no arg): real in [0, 1)."""
        return self._next_word() / 4294967296.0


# Hardcoded period-sampling polynomial coefficients (DERTY2.PAS lines 1401-1412),
# fit once (offline, via itgtr) against a catalog of real visual-binary orbits.
XFK = [
    -2.20038592302443e+0000,
     3.95320923444222e-0001,
    -1.40651047392788e-0002,
     2.97817055140045e-0004,
    -3.75276309676466e-0006,
     2.96148505471758e-0008,
    -1.51821284077818e-0010,
     5.13531862870052e-0013,
    -1.13748354489786e-0015,
     1.58736519209883e-0018,
    -1.26633323012809e-0021,
     4.40435754835436e-0025,
]


def get_ran_per2(rng):
    """Pascal GetRanPer2: rprd = random*500, period = exp(poly_11(rprd))."""
    rprd = rng.random() * 500.0
    nev = 0.0
    mn  = 1.0
    for coeff in XFK:
        nev += mn * coeff
        mn  *= rprd
    return math.exp(nev)


def get_ran_per(rng):
    """Pascal GetRanPer: GetRanPer2 with per-star rejection loops omitted
    (those only trigger for specific catalogued star designations)."""
    return get_ran_per2(rng)


def pascal_trial_search(rng, epochs, rho_obs, theta_obs, weights, excluded,
                        n_trials=2000):
    """One odflM4d random-trial sweep: draw n_trials candidates, score each
    with fit_thiele_innes (== GetPrmElOsn), return ascending-sorted by oev."""
    candidates = []
    for _ in range(n_trials):
        ev = rng.random()
        if ev < 0.99:
            prd = get_ran_per(rng)
            th  = rng.random() * prd - prd * 0.5
            res = fit_thiele_innes(epochs, rho_obs, theta_obs, weights, excluded, ev, th, prd)
            oev = res[4] if res is not None else 1e38
        else:
            prd = th = 0.0
            res = None
            oev = 1e5
        candidates.append((oev, ev, th, prd, res))

    candidates.sort(key=lambda c: c[0])   # stable ascending, matches Pascal sort7
    return candidates


def pascal_fit_pass(rng, epochs, rho_obs, theta_obs, weights, excluded,
                    n_trials=2000, n_refine=10):
    """One full odflM4d pass: 2000 trials -> top 10 refined -> best kept."""
    candidates = pascal_trial_search(rng, epochs, rho_obs, theta_obs, weights,
                                     excluded, n_trials)

    best_rms = 1e38;  best_params = None;  best_ti = None
    for oev, ev, th, prd, res in candidates[:n_refine]:
        if res is None:
            continue
        A, B, F, G, _ = res
        params, ti, rms = local_refine(
            [ev, th, prd, A, B, F, G], epochs, rho_obs, theta_obs, weights,
            excluded, max_nfev=20000
        )
        if ti is not None and rms < best_rms:
            best_rms = rms;  best_params = params;  best_ti = ti

    return best_params, best_ti, best_rms


# ---------------------------------------------------------------------------
# Full fitting pipeline
# ---------------------------------------------------------------------------

def fit_orbit(epochs_abs, rho_obs, theta_obs, verbose=True):
    """
    Complete orbit fitting.

    verbose — print progress messages (global search / outlier rounds) to stdout.

    Returns (ti, excluded, t0, epochs_rel).
    """
    def vprint(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    n  = len(epochs_abs)
    t0 = sum(epochs_abs) / n
    epochs = [e - t0 for e in epochs_abs]    # centred epochs

    excluded = [False] * n
    weights  = [1.0]   * n                   # uniform (no vesa3.res)

    vprint(f"\n  N_obs = {n},  t0 = {t0:.6f}", flush=True)

    best_rms    = 1e38
    best_params = None
    best_ti     = None

    # ------------------------------------------------------------------
    # Stage 0/1: Pascal-faithful global search (odflM4d, called twice in
    # TILIN132.PAS's main over one continuous RNG stream — the second
    # pass's winner is what the real program writes out).
    # ------------------------------------------------------------------
    rng = PascalRandom(seed=0)
    for pass_num in (1, 2):
        vprint(f"  Global random search (pass {pass_num}/2, 2000 trials) ...", flush=True)
        params, ti, rms = pascal_fit_pass(
            rng, epochs, rho_obs, theta_obs, weights, excluded
        )
        if ti is not None:
            vprint(f"    Pass {pass_num}: e={ti[1]:.4f}  P={ti[3]:.1f}  RMS={rms:.5f}",
                  flush=True)
            best_rms, best_params, best_ti = rms, params, ti

    if best_ti is None:
        raise RuntimeError("Global search found no valid orbit.")

    vprint(f"  Best solution so far: RMS = {best_rms:.5f}", flush=True)

    # ------------------------------------------------------------------
    # Stage 2: iterative outlier rejection
    # ------------------------------------------------------------------
    ti     = best_ti
    params = best_params

    for rnd in range(30):
        worst_k, worst_sigma = find_worst_outlier(ti, epochs, rho_obs, theta_obs, excluded)
        if worst_k < 0 or worst_sigma <= OUTLIER_SIGMA:
            break
        excluded[worst_k] = True
        n_used = sum(1 for ex in excluded if not ex)
        vprint(f"  Outlier round {rnd+1}: excluded obs #{worst_k} "
              f"(σ={worst_sigma:.1f}),  n_used={n_used}", flush=True)
        if n_used < 5:
            break

        params, ti, rms = local_refine(
            params, epochs, rho_obs, theta_obs, weights, excluded,
            max_nfev=20000
        )
        if ti is None:
            break
        vprint(f"    RMS = {rms:.5f}", flush=True)

    # ------------------------------------------------------------------
    # Stage 3: final polish
    # ------------------------------------------------------------------
    vprint("  Final polish ...", flush=True)
    params, ti, rms = local_refine(
        params, epochs, rho_obs, theta_obs, weights, excluded,
        max_nfev=200000, tol=1e-13
    )
    if ti is None:
        ti  = best_ti
        rms = best_rms
    else:
        vprint(f"  Final RMS = {rms:.6f}", flush=True)

    return ti, excluded, t0, epochs


# ---------------------------------------------------------------------------
# Read input.txt
# ---------------------------------------------------------------------------

def read_input(filename='input.txt'):
    epochs = [];  theta = [];  rho = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    epochs.append(float(parts[0]))
                    theta.append(float(parts[1]))
                    rho.append(float(parts[2]))
                except ValueError:
                    pass
    return epochs, theta, rho


# ---------------------------------------------------------------------------
# Radial velocity from Campbell elements  (Pascal: GetPvdFrmEl)
# ---------------------------------------------------------------------------

def compute_vr(abp, e, omm, prd, th, ai):
    """
    Line-of-sight velocity in km/s from Campbell elements, using finite
    differences on the z-coordinate of the 3-D position (Pascal GetPvdFrmEl).

    z(t) = abp * sin(ai) * (sin(omm)*(cos(E)-e) + cos(omm)*sqrt(1-e²)*sin(E))
    vr   = dz/dt / YEAR_TO_AU        [arcsec/yr treated as AU/yr at pi=1 → km/s]
    """
    sin_ai = math.sin(ai)
    som    = math.sin(omm)
    com    = math.cos(omm)
    e2     = math.sqrt(max(0.0, 1.0 - e * e))

    def z_at(t):
        M  = (t - th) * 2.0 * math.pi / prd
        E  = solve_kepler(M, e)
        xi = math.cos(E) - e
        yi = e2 * math.sin(E)
        return abp * sin_ai * (som * xi + com * yi)

    eps = 1e-8
    zt  = (z_at(+eps) - z_at(-eps)) / (2.0 * eps)
    return zt / YEAR_TO_AU


def _fnev(rho0, theta0, rho_dot, theta_dot, kriv, vr, mass,
          epochs_rel, rho_obs, theta_obs, excluded):
    """Pascal fnev: RMS residuals for a phase-space state (used by GetZnVr)."""
    if kriv < 0 or kriv > 3000 or mass < 0 or mass > 3000:
        return 1e38
    ti_t = tpvd(rho0, theta0, rho_dot, theta_dot, kriv, vr, mass, 0.0, 1.0)
    if ti_t is None:
        return 1e38
    if abs(ti_t[1] - 1.0) < 1e-9 or ti_t[1] > 10 or ti_t[3] < 1.0:
        return 1e38
    s = 0.0;  cnt = 0
    for k in range(len(epochs_rel)):
        if excluded[k]:
            continue
        try:
            rho_c, theta_c = rtti(ti_t, epochs_rel[k])
        except Exception:
            continue
        drho   = rho_obs[k] - rho_c
        dtheta = norm_deg(theta_obs[k] - theta_c)
        dte    = dtheta * math.pi * rho_obs[k] / 180.0
        s += drho * drho + dte * dte
        cnt += 1
    return math.sqrt(s / (2.0 * cnt)) if cnt > 0 else 1e38


def get_zn_vr(rho0, theta0, rho_dot, theta_dot, kriv, vr, mass,
              epochs_rel, rho_obs, theta_obs, excluded):
    """Pascal GetZnVr: pick the sign of vr that gives lower sky-plane RMS."""
    f0 = _fnev(rho0, theta0, rho_dot, theta_dot, kriv,  vr, mass,
               epochs_rel, rho_obs, theta_obs, excluded)
    f1 = _fnev(rho0, theta0, rho_dot, theta_dot, kriv, -vr, mass,
               epochs_rel, rho_obs, theta_obs, excluded)
    return vr if f0 <= f1 else -vr


# ---------------------------------------------------------------------------
# Derived quantities + O-C table from a fitted orbit — shared by write_output
# (file-based CLI) and the DataFrame-based API in tilin_module.py.
# ---------------------------------------------------------------------------

def compute_fit_summary(ti, excluded, t0, epochs_rel, epochs_abs,
                        rho_obs, theta_obs, weights):
    """
    Compute Campbell elements, derived quantities, and the O-C residual table
    for a fitted orbit state vector ti.

    Returns a dict with keys:
        t0, rho0, theta0, rho_dot, theta_dot, curvature, rv, mass_pi1,
        a, e, T, P, omega, Node, i, rms, rms_all, n_obs, n_excluded,
        rows  — list of per-observation dicts (epoch, rho_obs, theta_obs,
                 rho_c, theta_c, drho, dtheta, x, y, excluded), one per
                 input observation, in input order.
    """
    n    = len(epochs_abs)
    abp  = ti[0];  e = ti[1];  th = ti[2];  prd = ti[3]
    omm  = ti[8];  omb = ti[9];  ai = ti[10]

    # Inclination as in Pascal: tbk = 180 - ai*180/pi
    tbk = 180.0 - ai * 180.0 / math.pi
    if tbk < 0.0:
        tbk += 360.0

    # Initial-condition params at t=0 (relative to t0)
    # Numerically compute ρ₀, θ₀, dρ/dt, dθ/dt, curvature
    eps_t = 1e-4
    try:
        rho0_fit,   theta0_fit   = rtti(ti, 0.0)
        rho_m, theta_m = rtti(ti, -eps_t)
        rho_p, theta_p = rtti(ti,  eps_t)
        rho_dot   = (rho_p  - rho_m)  / (2.0 * eps_t)
        theta_dot = norm_deg(theta_p - theta_m) / (2.0 * eps_t)
    except Exception:
        rho0_fit = abp;  theta0_fit = 0.0;  rho_dot = 0.0;  theta_dot = 0.0

    # Curvature at t=0 via second derivative of sky position
    try:
        eps2 = 1e-3
        r_a, t_a = rtti(ti, -eps2)
        r_b, t_b = rtti(ti,  eps2)
        r_c, t_c = rtti(ti,   0.0)
        xa_  = r_a * math.sin(t_a * DEG);  ya_ = r_a * math.cos(t_a * DEG)
        xb_  = r_b * math.sin(t_b * DEG);  yb_ = r_b * math.cos(t_b * DEG)
        xc_  = r_c * math.sin(t_c * DEG);  yc_ = r_c * math.cos(t_c * DEG)
        vx0  = (xb_ - xa_) / (2*eps2);      vy0 = (yb_ - ya_) / (2*eps2)
        ax0  = (xa_ - 2*xc_ + xb_) / eps2**2
        ay0  = (ya_ - 2*yc_ + yb_) / eps2**2
        mu   = math.sqrt(vx0**2 + vy0**2)
        cross = vx0*ay0 - vy0*ax0
        kriv  = mu**2 / abs(cross) if abs(cross) > 1e-30 else 0.0
    except Exception:
        kriv = 0.0

    # Mass (a³/P²) for parallax=1 — same as Pascal tmassa when mbcir=false
    mass_pi1 = abp**3 / prd**2 if prd > 0.0 else 0.0

    # RMS — used observations (kvr=1 since all weights=1 and sumves=nt2)
    n_used   = sum(1 for ex in excluded if not ex)
    sum_used = 0.0
    sum_all  = 0.0
    rows     = []
    for k in range(n):
        try:
            rho_c, theta_c = rtti(ti, epochs_rel[k])
        except Exception:
            rho_c = 0.0;  theta_c = 0.0
        drho   = rho_obs[k] - rho_c
        dtheta = norm_deg(theta_obs[k] - theta_c)
        dte    = dtheta * math.pi * rho_obs[k] / 180.0
        w      = weights[k]
        contrib = w*w * (drho*drho + dte*dte)
        sum_all += contrib
        if not excluded[k]:
            sum_used += contrib
        x_obs = rho_obs[k] * math.sin(theta_obs[k] * DEG)
        y_obs = rho_obs[k] * math.cos(theta_obs[k] * DEG)
        rows.append({
            'epoch':     epochs_abs[k],
            'rho_obs':   rho_obs[k],
            'theta_obs': theta_obs[k],
            'rho_c':     rho_c,
            'theta_c':   theta_c,
            'drho':      drho,
            'dtheta':    dtheta,
            'x':         x_obs,
            'y':         y_obs,
            'excluded':  excluded[k],
        })

    rms_used = math.sqrt(0.5 * sum_used / n_used) if n_used > 0 else 0.0
    rms_all  = math.sqrt(0.5 * sum_all  / n)      if n      > 0 else 0.0

    vr = compute_vr(abp, e, omm, prd, th, ai)
    vr = get_zn_vr(rho0_fit, theta0_fit, rho_dot, theta_dot, kriv, vr,
                   mass_pi1, epochs_rel, rho_obs, theta_obs, excluded)

    return {
        't0':         t0,
        'rho0':       rho0_fit,
        'theta0':     theta0_fit,
        'rho_dot':    rho_dot,
        'theta_dot':  theta_dot,
        'curvature':  kriv,
        'rv':         vr,
        'mass_pi1':   mass_pi1,
        'a':          abp,
        'e':          e,
        'T':          th + t0,
        'P':          prd,
        'omega':      omm * 180.0 / math.pi,
        'Node':       omb * 180.0 / math.pi,
        'i':          tbk,
        'rms':        rms_used,
        'rms_all':    rms_all,
        'n_obs':      n,
        'n_excluded': n - n_used,
        'rows':       rows,
    }


# ---------------------------------------------------------------------------
# Write output.txt  (matches vydres4 format exactly)
# ---------------------------------------------------------------------------

def write_output(ti, excluded, t0, epochs_rel, epochs_abs,
                 rho_obs, theta_obs, weights, filename='output.txt'):
    """Write results in the format of the original Pascal program."""
    s = compute_fit_summary(ti, excluded, t0, epochs_rel, epochs_abs,
                            rho_obs, theta_obs, weights)

    with open(filename, 'w') as f:
        f.write(' \n')
        f.write(f"t0            ={s['t0']:22.15f}\n")
        f.write(f"rho0          ={s['rho0']:22.15f}\n")
        f.write(f"theta0        ={s['theta0']:22.15f}\n")
        f.write(f"rho t         ={s['rho_dot']:22.15f}\n")
        f.write(f"theta t       ={s['theta_dot']:22.15f}\n")
        f.write(f"rad.curvature ={s['curvature']:22.15f}\n")
        f.write(f"R V[\"/year]   ={s['rv']:22.15f}\n")
        f.write(f"mass(pi=1)    ={s['mass_pi1']:22.15f}\n")
        f.write('\n')
        f.write(f"a    ={s['a']:22.15f}\n")
        f.write(f"e    ={s['e']:22.15f}\n")
        f.write(f"T    ={s['T']:22.15f}\n")
        f.write(f"P    ={s['P']:22.15f}\n")
        f.write(f"omega={s['omega']:22.15f}\n")
        f.write(f"Node ={s['Node']:22.15f}\n")
        f.write(f"i    ={s['i']:22.15f}\n")
        f.write('\n')
        f.write(f"RMS={s['rms']:22.15f}\n")
        f.write('\n')
        f.write('RMS\n')
        f.write(f"all obs.={s['rms_all']:22.15f}\n")
        f.write('\n')
        f.write('    Ep       Rho(obs)      Theta(obs)     Rho(c)     Theta(c)       dRho        dTheta          X            Y \n')

        for row in s['rows']:
            if row['excluded']:
                continue
            f.write(f"{row['epoch']:10.5f} {row['rho_obs']:12.7f} {row['theta_obs']:12.7f}"
                    f" {row['rho_c']:12.7f} {row['theta_c']:12.7f}"
                    f" {row['drho']:12.7f} {row['dtheta']:12.7f} {row['x']:12.7f} {row['y']:12.7f}\n")

        f.write('\n')
        f.write('excluded:\n')

        for row in s['rows']:
            if not row['excluded']:
                continue
            f.write(f"{row['epoch']:10.5f} {row['rho_obs']:12.7f} {row['theta_obs']:12.7f}"
                    f" {row['rho_c']:12.7f} {row['theta_c']:12.7f}"
                    f" {row['drho']:12.7f} {row['dtheta']:12.7f} {row['x']:12.7f} {row['y']:12.7f}\n")

        f.write('\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("TILIN132 — Visual binary orbit fitting  (Python translation)")
    print("=" * 62)

    epochs_abs, theta_obs, rho_obs = read_input('input.txt')
    if len(epochs_abs) < 3:
        print("ERROR: need at least 3 observations.", file=sys.stderr)
        sys.exit(1)

    n = len(epochs_abs)
    print(f"  Observations : {n}")
    print(f"  Epoch range  : {min(epochs_abs):.2f} – {max(epochs_abs):.2f}")

    ti, excluded, t0, epochs_rel = fit_orbit(epochs_abs, rho_obs, theta_obs)
    campbell = campbell_from_ti(ti, t0)

    n_excl = sum(excluded)
    rms    = rms_residual(ti, epochs_rel, rho_obs, theta_obs,
                          [1.0]*n, excluded)

    print()
    print("=" * 62)
    print(f"  a     = {campbell['a']:.6f} arcsec")
    print(f"  e     = {campbell['e']:.6f}")
    print(f"  T     = {campbell['T']:.3f} yr")
    print(f"  P     = {campbell['P']:.3f} yr")
    print(f"  ω     = {campbell['omega']:.3f} deg")
    print(f"  Ω     = {campbell['Node']:.3f} deg")
    print(f"  i     = {campbell['i']:.3f} deg")
    print(f"  RMS   = {rms:.6f} arcsec")
    print(f"  Excl. = {n_excl}/{n}")
    print("=" * 62)

    write_output(
        ti, excluded, t0, epochs_rel, epochs_abs,
        rho_obs, theta_obs, [1.0]*n,
        filename='output.txt'
    )
    print("Results written to output.txt")


if __name__ == '__main__':
    main()
