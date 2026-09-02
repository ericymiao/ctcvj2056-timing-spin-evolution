"""One-shot edit of notebook 08 to unify the IP propeller/accretion regime.

Changes:
  1. Constants cell: remove module-level xi = 0.5; keep eta = 1.0.
  2. Forward-models markdown: update propeller description.
  3. Forward-models code: thread xi through pdot_ip / kT_ip / LX_ip.
  4. Priors/likelihood markdown: document 4D IP prior and branch-dependent
     likelihood.
  5. Priors/likelihood code: 4D prior_transform_ip, unified log_likelihood_ip
     with accretion/propeller branches.
  6. IP nested sampling: ndim 3 -> 4.
  7. Summary cell: include xi in IP labels.
  8. NEW cell (regime-partitioned summary) inserted after the summary cell.
  9. IP corner plots: marginalise over xi (3D plot retained).
 10. Predicted Pdot / PPC cells: pass xi to pdot_ip / kT_ip / LX_ip.
 11. Old-Pdot comparison: rebuild log_likelihood_ip_old and sampler to 4D.
 12. Comparison table & corner overlay: 4D-aware.

Run this script ONCE to update the notebook in place.
"""
import json
import pathlib
import shutil

nb_path = pathlib.Path(__file__).with_name("08_bayesian_spectral_constraints.ipynb")
backup = nb_path.with_suffix(".ipynb.bak")
if not backup.exists():
    shutil.copy2(nb_path, backup)

nb = json.loads(nb_path.read_text())
cells = nb["cells"]


def set_source(idx, text):
    cells[idx]["source"] = text.splitlines(keepends=True)
    if cells[idx]["cell_type"] == "code":
        cells[idx]["outputs"] = []
        cells[idx]["execution_count"] = None


# ---------------------------------------------------------------------------
# Cell 3 — constants (remove module-level xi)
# ---------------------------------------------------------------------------
set_source(3, '''# --- Physical constants (CGS) ---
G = 6.674e-8       # cm^3 g^-1 s^-2
Msun = 1.989e33    # g
Rsun = 6.957e10    # cm
m_H = 1.6726e-24   # proton mass [g]
keV_to_erg = 1.602e-9  # keV -> erg

# --- Measured timing data (1-sigma) ---
P = 29.609685843       # spin period [s]
Pdot_obs = 0.0         # observed Pdot [s/s]
sigma_Pdot = 1.7e-15   # 1-sigma uncertainty on Pdot [s/s]

# --- Spectral observables (Salcedo+2024) ---
# X-ray luminosity: F_X(unabs, 0.3-12 keV) = 2.2e-12 erg/s/cm^2, d = 261.6 +/- 7.4 pc
d_pc = 261.6                           # distance [pc]
sigma_d_pc = 7.4                       # 1-sigma distance error [pc]
d_cm = d_pc * 3.086e18                 # distance [cm]
F_X = 2.2e-12                          # unabsorbed flux [erg/s/cm^2]
sigma_FX_frac = 0.05                   # ~5% flux uncertainty (spectral fit)
L_X_obs = 4 * np.pi * d_cm**2 * F_X   # X-ray luminosity [erg/s], 0.3-12 keV band
# Propagate distance + flux uncertainties
sigma_LX = L_X_obs * np.sqrt((2 * sigma_d_pc / d_pc)**2 + sigma_FX_frac**2)

# X-ray band for the bolometric -> band-limited correction applied to L_X_pred
LX_BAND_LO = 0.3   # keV
LX_BAND_HI = 12.0  # keV

# kT_max from MKCFLOW (Salcedo+2024, Table 6)
# 90% CI: 16.3 +1.9/-1.8 keV -> 1-sigma ~ 1.9/1.645 ~ 1.15 keV
kT_obs = 16.3        # keV
sigma_kT = 1.15       # 1-sigma [keV]

# Mean molecular weight for fully ionized solar-composition plasma
mu_mol = 0.617
# Empirical correction factor (Yu+2018, calibrated on 11 nmCVs)
alpha_Yu = 0.646

# --- Load WD mass-radius-inertia table ---
mr_data = np.loadtxt('../wd_mass_radius/mass_radius.txt', delimiter=',', skiprows=1)
M_table = mr_data[:, 0]   # M_sun
R_table = mr_data[:, 1]   # R_sun
I_table = mr_data[:, 2]   # g cm^2

# Interpolators: input in solar masses, output in CGS
get_R = interp1d(M_table, R_table * Rsun, kind='cubic')   # returns cm
get_I = interp1d(M_table, I_table, kind='cubic')           # returns g cm^2

M_min_table = M_table.min()
M_max_table = M_table.max()

# --- Magnetosphere geometry ---
# xi:  magnetospheric-radius prefactor, r_m = xi * r_Alfven (Ghosh & Lamb 1979).
#      Conventional disk-accretion value is ~0.5, literature range ~0.3-1.
#      Since r_m ~ B^(4/7), the inferred B scales as xi^(-7/4); a factor-2
#      uncertainty in xi maps to ~3.4x on B. We marginalise over xi as a
#      sampled nuisance parameter so the B posterior reflects this systematic
#      directly (see prior section below).
# eta: accretion-efficiency factor in the propeller regime (Bhattacharyya &
#      Chakrabarty 2017, their Eq 10). eta = 1 assumes full coupling between
#      ejected material and the magnetosphere. Held fixed; the posterior
#      fraction in the propeller branch depends only weakly on eta.
eta = 1.0

print(f"Table range: {M_min_table:.3f} - {M_max_table:.3f} Msun")
print(f"L_X_obs = {L_X_obs:.2e} erg/s  (sigma = {sigma_LX:.2e})")
print(f"kT_obs  = {kT_obs:.1f} +/- {sigma_kT:.2f} keV (1-sigma)")
''')


# ---------------------------------------------------------------------------
# Cell 4 — forward-models markdown (update propeller description)
# ---------------------------------------------------------------------------
set_source(4, '''## Forward Models

**Pdot** models implement Bhattacharyya & Chakrabarty (2017, Eqs. 9-10). The IP torque uses a propeller-efficiency factor $\\eta = 1$.

**Spectral forward models:**

*nmCV:*
- $kT_{\\text{pred}} = \\frac{3}{16} \\alpha \\mu_{\\text{mol}} m_H \\frac{GM}{R}$ (Yu+2018, $\\alpha = 0.646$)
- $L_{X,\\text{pred}} = f_{\\rm band}(kT_{\\rm pred}) \\cdot \\frac{GM\\dot{M}}{R}$

*IP (accretion branch only, $r_m \\le r_{\\rm co}$):*
- $kT_{\\text{pred}} = \\frac{3}{8} \\mu_{\\text{mol}} m_H \\cdot 2GM\\left(\\frac{1}{R} - \\frac{1}{r_m}\\right)$ (compact shock approximation)
- $L_{X,\\text{pred}} = f_{\\rm band}(kT_{\\rm pred}) \\cdot GM\\dot{M}\\left(\\frac{1}{R} - \\frac{1}{r_m}\\right)$

*IP (propeller branch, $r_m > r_{\\rm co}$):* no standing accretion shock on the WD — $kT_{\\rm ip}$ / $L_{X,\\rm ip}$ return NaN. The likelihood (below) drops the $kT$ term in this branch and treats the observed $L_X$ as an upper limit on any residual accretion luminosity $\\sim GM\\dot M/r_m$.

**Band correction.** Theoretical accretion luminosities are bolometric, whereas $L_{X,\\rm obs}$ is the unabsorbed flux in 0.3-12 keV. For single-temperature thermal bremsstrahlung with Gaunt factor $g_{\\rm ff} \\approx 1$,
$$f_{\\rm band}(kT) = \\exp(-E_{\\rm lo}/kT) - \\exp(-E_{\\rm hi}/kT),$$
evaluated at the predicted $kT$. At $kT \\sim 16$ keV this gives $f_{\\rm band} \\approx 0.5$.
''')


# ---------------------------------------------------------------------------
# Cell 5 — forward models code (thread xi as argument)
# ---------------------------------------------------------------------------
set_source(5, '''# =====================================================
# Pdot forward models (Bhattacharyya & Chakrabarty 2017)
# =====================================================

def pdot_nmcv(M_wd_solar, log10_Mdot):
    """Predicted Pdot for the nmCV case (disk extends to WD surface, spin-up only)."""
    M = M_wd_solar * Msun
    Mdot = 10**log10_Mdot
    R = get_R(M_wd_solar)
    I = get_I(M_wd_solar)
    tau = Mdot * np.sqrt(G * M * R)
    return -P**2 * tau / (2 * np.pi * I)


def pdot_ip(M_wd_solar, log10_Mdot, log10_mu, xi):
    """Predicted Pdot for IP case (piecewise disk-magnetosphere torque).

    Accretion phase (r_m <= r_co): Eq 9 of Bhattacharyya & Chakrabarty (2017).
    Propeller phase (r_m > r_co):  Eq 10 with efficiency `eta` (module-level).
    Returns NaN only if r_m <= R_WD (magnetosphere crushed to WD surface).
    """
    M = M_wd_solar * Msun
    Mdot = 10**log10_Mdot
    mu = 10**log10_mu
    R = get_R(M_wd_solar)
    I = get_I(M_wd_solar)

    r_m = xi * (mu**4 / (2 * G * M * Mdot**2))**(1.0/7.0)
    r_co = (G * M * P**2 / (4 * np.pi**2))**(1.0/3.0)

    if r_m <= R:
        return np.nan

    if r_m <= r_co:
        x = r_m / r_co
        tau_acc = Mdot * np.sqrt(G * M * r_m)
        tau_mag = mu**2 / (9 * r_m**3) * (2 * x**3 - 6 * x**1.5 + 3)
    else:
        tau_acc = -eta * Mdot * np.sqrt(G * M * r_m)
        tau_mag = -mu**2 / (9 * r_m**3) * (3 - 2 * (r_co / r_m)**1.5)

    return -P**2 * (tau_acc + tau_mag) / (2 * np.pi * I)


# =====================================================
# Band-correction factor for the bolometric -> 0.3-12 keV conversion
# =====================================================

def f_band(kT_keV):
    """Fraction of bolometric thermal-bremsstrahlung luminosity in LX_BAND_LO..LX_BAND_HI.

    Single-temperature bremsstrahlung with Gaunt factor ~1:
        dL/dE ~ exp(-E/kT)  ->  L(E1,E2)/L_bol = exp(-E1/kT) - exp(-E2/kT).
    """
    return np.exp(-LX_BAND_LO / kT_keV) - np.exp(-LX_BAND_HI / kT_keV)


# =====================================================
# Spectral forward models
# =====================================================

def kT_nmcv(M_wd_solar):
    """Predicted maximum shock temperature for nmCV [keV] (Yu+2018 relation)."""
    M = M_wd_solar * Msun
    R = get_R(M_wd_solar)
    kT_erg = (3.0 / 16.0) * alpha_Yu * mu_mol * m_H * G * M / R
    return kT_erg / keV_to_erg


def LX_nmcv(M_wd_solar, log10_Mdot):
    """Predicted 0.3-12 keV X-ray luminosity for nmCV [erg/s]."""
    M = M_wd_solar * Msun
    R = get_R(M_wd_solar)
    Mdot = 10**log10_Mdot
    L_bol = G * M * Mdot / R
    return f_band(kT_nmcv(M_wd_solar)) * L_bol


def kT_ip(M_wd_solar, log10_Mdot, log10_mu, xi):
    """Predicted shock temperature for IP [keV], accretion branch only.

    Returns NaN if r_m <= R_WD (unphysical) or r_m > r_co (propeller: no
    standing accretion shock on the WD).
    """
    M = M_wd_solar * Msun
    R = get_R(M_wd_solar)
    Mdot = 10**log10_Mdot
    mu = 10**log10_mu
    r_m = xi * (mu**4 / (2 * G * M * Mdot**2))**(1.0/7.0)
    r_co = (G * M * P**2 / (4 * np.pi**2))**(1.0/3.0)
    if r_m <= R or r_m > r_co:
        return np.nan
    vff2 = 2 * G * M * (1.0 / R - 1.0 / r_m)
    kT_erg = (3.0 / 8.0) * mu_mol * m_H * vff2
    return kT_erg / keV_to_erg


def LX_ip(M_wd_solar, log10_Mdot, log10_mu, xi):
    """Predicted 0.3-12 keV X-ray luminosity for IP [erg/s], accretion branch."""
    M = M_wd_solar * Msun
    R = get_R(M_wd_solar)
    Mdot = 10**log10_Mdot
    mu = 10**log10_mu
    r_m = xi * (mu**4 / (2 * G * M * Mdot**2))**(1.0/7.0)
    r_co = (G * M * P**2 / (4 * np.pi**2))**(1.0/3.0)
    if r_m <= R or r_m > r_co:
        return np.nan
    vff2 = 2 * G * M * (1.0 / R - 1.0 / r_m)
    kT = (3.0 / 8.0) * mu_mol * m_H * vff2 / keV_to_erg
    L_bol = G * M * Mdot * (1.0 / R - 1.0 / r_m)
    return f_band(kT) * L_bol


# --- Sanity checks ---
print("=== Pdot checks ===")
print(f"nmCV(0.935, 15): Pdot = {pdot_nmcv(0.935, 15):.3e}")
print(f"IP(0.815, 15, 30.3, xi=0.5): Pdot = {pdot_ip(0.815, 15, 30.3, 0.5):.3e}")
# AE-Aqr-like propeller point should give strong negative Pdot, not NaN:
print(f"IP(0.8, 14, 33, xi=0.5) [propeller]: Pdot = {pdot_ip(0.8, 14, 33, 0.5):.3e}")

print("\\n=== f_band checks ===")
print(f"f_band(16.3 keV): {f_band(16.3):.3f}   (0.3-12 keV fraction at kT_obs)")
print(f"f_band(10 keV):   {f_band(10.0):.3f}")
print(f"f_band(30 keV):   {f_band(30.0):.3f}")

print("\\n=== kT checks ===")
print(f"nmCV kT(0.9): {kT_nmcv(0.9):.1f} keV  (expect ~15-20 keV)")
print(f"IP kT(0.8, 14.5, 30.0, 0.5): {kT_ip(0.8, 14.5, 30.0, 0.5):.1f} keV  (accretion)")
kT_prop = kT_ip(0.8, 14.5, 31.0, 0.5)
print(f"IP kT(0.8, 14.5, 31.0, 0.5): {kT_prop}   (propeller -> NaN expected)")

print("\\n=== L_X checks (band-corrected) ===")
print(f"nmCV L_X(0.9, 14.5): {LX_nmcv(0.9, 14.5):.2e} erg/s  (expect ~5e30-1e31)")
print(f"IP L_X(0.8, 14.5, 30.0, 0.5): {LX_ip(0.8, 14.5, 30.0, 0.5):.2e} erg/s")
''')


# ---------------------------------------------------------------------------
# Cell 6 — priors/likelihood markdown
# ---------------------------------------------------------------------------
set_source(6, '''## Priors and Likelihood

**Priors (nmCV, 2D):**
- $M_{\\rm WD}$: $\\mathcal{N}(0.83, 0.20)$ (Pala+2022, volume-limited CV sample), clipped to WD table range
- $\\log_{10}\\dot{M}$: $\\mathcal{U}(13, 17)$ (flat in log-space)

**Priors (IP, 4D):**
- $M_{\\rm WD}$: $\\mathcal{N}(0.83, 0.20)$ clipped to table range
- $\\log_{10}\\dot{M}$: $\\mathcal{U}(13, 17)$
- $\\log_{10}\\mu$: $\\mathcal{U}(29, 35)$ (magnetic dipole moment, G$\\,$cm$^3$)
- $\\xi$: log-uniform on $[0.3, 1.0]$ (magnetospheric-radius prefactor; Ghosh & Lamb 1979 range)

Sampling $\\xi$ marginalises the factor-$\\xi^{-7/4}$ systematic in $B$ into the posterior rather than fixing it at a conventional value.

**Likelihood (branch-dependent for IP).** All three observables enter as independent Gaussians where applicable:
$$\\chi^2_{\\dot P} = \\left(\\frac{\\dot P_{\\rm obs} - \\dot P_{\\rm pred}}{\\sigma_{\\dot P}}\\right)^2,\\quad \\chi^2_{kT} = \\left(\\frac{kT_{\\rm obs} - kT_{\\rm pred}}{\\sigma_{kT}}\\right)^2,\\quad \\chi^2_{L_X} = \\left(\\frac{L_{X,\\rm obs} - L_{X,\\rm pred}}{\\sigma_{L_X}}\\right)^2$$

For the nmCV model: $\\ln\\mathcal{L} = -\\tfrac12(\\chi^2_{\\dot P} + \\chi^2_{kT} + \\chi^2_{L_X})$.

For the IP model:
- **Accretion branch** ($r_m \\le r_{\\rm co}$): same as nmCV — all three $\\chi^2$ terms contribute.
- **Propeller branch** ($r_m > r_{\\rm co}$): $\\chi^2_{\\dot P}$ always included; $\\chi^2_{kT}$ dropped (no accretion shock on WD); $L_X$ treated as a 2-$\\sigma$ upper limit on any residual accretion luminosity $L_{X,\\rm resid} = \\eta GM\\dot M/r_m$ — a half-Gaussian penalty fires only when $L_{X,\\rm resid} > L_{X,\\rm obs} + 2\\sigma_{L_X}$.

This lets the sampler enter the propeller regime; whether it stays there is decided by the $\\dot P$ data, not by the likelihood code.
''')


# ---------------------------------------------------------------------------
# Cell 7 — priors/likelihood code
# ---------------------------------------------------------------------------
set_source(7, '''# =====================================================
# nmCV model: 2 parameters [M_wd, log10_Mdot]
# =====================================================

# Population mass prior (Pala+2022)
POP_M_MU, POP_M_SIG = 0.83, 0.20

# log10(Mdot) flat prior
LOGMDOT_LO, LOGMDOT_HI = 13.0, 17.0


def prior_transform_nmcv(u):
    """Map unit cube [0,1]^2 -> nmCV parameter space."""
    theta = np.empty(2)
    theta[0] = norm.ppf(u[0], loc=POP_M_MU, scale=POP_M_SIG)
    theta[0] = np.clip(theta[0], M_min_table + 0.001, M_max_table - 0.001)
    theta[1] = LOGMDOT_LO + u[1] * (LOGMDOT_HI - LOGMDOT_LO)
    return theta


def log_likelihood_nmcv(theta):
    """Log-likelihood for nmCV model with Pdot + L_X + kT."""
    M_wd, log_Mdot = theta
    if M_wd <= M_min_table or M_wd >= M_max_table:
        return -1e100

    Pdot_pred = pdot_nmcv(M_wd, log_Mdot)
    kT_pred = kT_nmcv(M_wd)
    LX_pred = LX_nmcv(M_wd, log_Mdot)

    chi2  = ((Pdot_obs - Pdot_pred) / sigma_Pdot)**2
    chi2 += ((kT_obs - kT_pred) / sigma_kT)**2
    chi2 += ((L_X_obs - LX_pred) / sigma_LX)**2
    return -0.5 * chi2


# =====================================================
# IP model: 4 parameters [M_wd, log10_Mdot, log10_mu, xi]
# =====================================================

LOGMU_LO, LOGMU_HI = 29.0, 35.0
XI_LO, XI_HI = 0.3, 1.0  # log-uniform prior range (Ghosh & Lamb 1979)


def prior_transform_ip(u):
    """Map unit cube [0,1]^4 -> IP parameter space."""
    theta = np.empty(4)
    theta[0] = norm.ppf(u[0], loc=POP_M_MU, scale=POP_M_SIG)
    theta[0] = np.clip(theta[0], M_min_table + 0.001, M_max_table - 0.001)
    theta[1] = LOGMDOT_LO + u[1] * (LOGMDOT_HI - LOGMDOT_LO)
    theta[2] = LOGMU_LO + u[2] * (LOGMU_HI - LOGMU_LO)
    # log-uniform on xi in [XI_LO, XI_HI]
    theta[3] = 10**(np.log10(XI_LO) + u[3] * (np.log10(XI_HI) - np.log10(XI_LO)))
    return theta


def log_likelihood_ip(theta):
    """Unified IP log-likelihood spanning accretion and propeller branches.

    Accretion branch (r_m <= r_co): Pdot + kT + L_X (all as Gaussians).
    Propeller  branch (r_m >  r_co): Pdot only; L_X as an upper limit on the
        residual accretion luminosity eta * G M Mdot / r_m (penalty fires
        only if that prediction exceeds L_X_obs + 2 sigma_LX).
    """
    M_wd, log_Mdot, log_mu, xi_s = theta
    if M_wd <= M_min_table or M_wd >= M_max_table:
        return -1e100

    # Predicted Pdot — valid in both branches (NaN only if r_m <= R_WD)
    Pdot_pred = pdot_ip(M_wd, log_Mdot, log_mu, xi_s)
    if np.isnan(Pdot_pred):
        return -1e100

    chi2 = ((Pdot_obs - Pdot_pred) / sigma_Pdot)**2

    # Decide which regime we are in from r_m / r_co
    M_cgs = M_wd * Msun
    Mdot = 10**log_Mdot
    mu = 10**log_mu
    r_m = xi_s * (mu**4 / (2 * G * M_cgs * Mdot**2))**(1.0/7.0)
    r_co = (G * M_cgs * P**2 / (4 * np.pi**2))**(1.0/3.0)

    if r_m <= r_co:
        # Accretion branch — standing shock on WD, all spectral terms apply
        kT_pred = kT_ip(M_wd, log_Mdot, log_mu, xi_s)
        LX_pred = LX_ip(M_wd, log_Mdot, log_mu, xi_s)
        if np.isnan(kT_pred) or kT_pred <= 0:
            return -1e100
        if np.isnan(LX_pred) or LX_pred <= 0:
            return -1e100
        chi2 += ((kT_obs - kT_pred) / sigma_kT)**2
        chi2 += ((L_X_obs - LX_pred) / sigma_LX)**2
    else:
        # Propeller branch — no accretion shock on WD.
        # Drop kT term. Use L_X as a 2-sigma upper limit on residual accretion.
        LX_resid = eta * G * M_cgs * Mdot / r_m
        LX_cap = L_X_obs + 2 * sigma_LX
        if LX_resid > LX_cap:
            chi2 += ((LX_resid - LX_cap) / sigma_LX)**2

    return -0.5 * chi2


# --- Quick check that the propeller regime is reachable ---
test_acc  = log_likelihood_ip([0.8, 14.5, 30.0, 0.5])
test_prop = log_likelihood_ip([0.8, 14.0, 33.0, 0.5])
print(f"log L(accretion point):  {test_acc:.3e}   (should be finite, typically O(1-100))")
print(f"log L(propeller  point): {test_prop:.3e}  (should be finite, large negative from Pdot chi2)")
''')


# ---------------------------------------------------------------------------
# Cell 9 — IP sampler: ndim 3 -> 4
# ---------------------------------------------------------------------------
set_source(9, '''print("Running IP nested sampling (4 params: M_WD, log Mdot, log mu, xi)...")
sampler_ip = dynesty.DynamicNestedSampler(
    log_likelihood_ip, prior_transform_ip, ndim=4,
    nlive=2000, bound='multi', sample='rwalk'
)
sampler_ip.run_nested(print_progress=True)
results_ip = sampler_ip.results
print(f"IP done: logZ = {results_ip['logz'][-1]:.2f} +/- {results_ip['logzerr'][-1]:.2f}")
''')


# ---------------------------------------------------------------------------
# Cell 12 — summary statistics (add xi)
# ---------------------------------------------------------------------------
set_source(12, '''# --- Extract equal-weight posterior samples ---

# IP
weights_ip = np.exp(results_ip['logwt'] - results_ip['logz'][-1])
weights_ip /= weights_ip.sum()
samples_ip_eq = dyfunc.resample_equal(results_ip['samples'], weights_ip)

# Derived: B = mu / R^3 (R depends on sampled M_WD)
R_samples_ip = np.array([get_R(m) for m in samples_ip_eq[:, 0]])
mu_samples_ip = 10**samples_ip_eq[:, 2]
B_samples_ip = mu_samples_ip / R_samples_ip**3
log10_B_samples_ip = np.log10(B_samples_ip)

# nmCV
weights_nmcv = np.exp(results_nmcv['logwt'] - results_nmcv['logz'][-1])
weights_nmcv /= weights_nmcv.sum()
samples_nmcv_eq = dyfunc.resample_equal(results_nmcv['samples'], weights_nmcv)


# --- Summary statistics ---
def summarize(samples, labels):
    for i, label in enumerate(labels):
        med = np.median(samples[:, i])
        lo, hi = np.percentile(samples[:, i], [16, 84])
        print(f'  {label}: {med:.4f}  (+{hi-med:.4f} / -{med-lo:.4f})')

print('=== IP posteriors ===')
summarize(samples_ip_eq,
          ['M_WD [Msun]', 'log10(Mdot) [g/s]', 'log10(mu) [G cm^3]', 'xi'])
med_B = np.median(log10_B_samples_ip)
lo_B, hi_B = np.percentile(log10_B_samples_ip, [16, 84])
print(f'  log10(B) [G]: {med_B:.2f}  (+{hi_B-med_B:.2f} / -{med_B-lo_B:.2f})')
print(f'  B [MG]: {10**med_B/1e6:.4f}')
print(f'  B [kG]: {10**med_B/1e3:.2f}')

print('\\n=== nmCV posteriors ===')
summarize(samples_nmcv_eq, ['M_WD [Msun]', 'log10(Mdot) [g/s]'])
''')


# ---------------------------------------------------------------------------
# Insert new markdown + code cell for regime-partitioned summary AFTER cell 12.
# After insertion the subsequent cells shift by +2.
# ---------------------------------------------------------------------------
regime_md = {
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "### Regime-partitioned posterior (IP)\n",
        "\n",
        "Split the IP posterior into **accretion** ($r_m \\le r_{\\rm co}$) and **propeller** ($r_m > r_{\\rm co}$) branches and report the fraction of posterior mass in each. This is the direct, data-driven answer to \u201cis CTCV J2056 an accreting IP or a propeller?\u201d\n"
    ],
}

regime_code = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": (
        "# --- Partition IP posterior into accretion (r_m <= r_co) and propeller (r_m > r_co) ---\n"
        "M_post      = samples_ip_eq[:, 0]\n"
        "logMdot_post = samples_ip_eq[:, 1]\n"
        "logmu_post   = samples_ip_eq[:, 2]\n"
        "xi_post      = samples_ip_eq[:, 3]\n"
        "\n"
        "M_cgs_post  = M_post * Msun\n"
        "Mdot_post   = 10**logMdot_post\n"
        "mu_post     = 10**logmu_post\n"
        "r_m_post  = xi_post * (mu_post**4 / (2 * G * M_cgs_post * Mdot_post**2))**(1.0/7.0)\n"
        "r_co_post = (G * M_cgs_post * P**2 / (4 * np.pi**2))**(1.0/3.0)\n"
        "\n"
        "acc_mask  = r_m_post <= r_co_post\n"
        "prop_mask = ~acc_mask\n"
        "f_acc  = acc_mask.mean()\n"
        "f_prop = prop_mask.mean()\n"
        "\n"
        "# Predicted Pdot at each posterior sample (both branches)\n"
        "Pdot_ip_samples = np.array([\n"
        "    pdot_ip(s[0], s[1], s[2], s[3]) for s in samples_ip_eq\n"
        "])\n"
        "\n"
        "print(f'Posterior fraction in accretion branch (r_m <= r_co): {f_acc:.3f}')\n"
        "print(f'Posterior fraction in propeller branch (r_m >  r_co): {f_prop:.3f}')\n"
        "\n"
        "def _branch_stats(mask, name):\n"
        "    n = int(mask.sum())\n"
        "    print(f'\\n=== {name} branch ({n} samples) ===')\n"
        "    if n == 0:\n"
        "        print('  (empty)')\n"
        "        return\n"
        "    for i, lbl in enumerate(['M_WD [Msun]', 'log10(Mdot) [g/s]',\n"
        "                             'log10(mu) [G cm^3]', 'xi']):\n"
        "        col = samples_ip_eq[mask, i]\n"
        "        med = np.median(col); lo, hi = np.percentile(col, [16, 84])\n"
        "        print(f'  {lbl:<22s} {med:.3f}  (+{hi-med:.3f}/-{med-lo:.3f})')\n"
        "    col = log10_B_samples_ip[mask]\n"
        "    med = np.median(col); lo, hi = np.percentile(col, [16, 84])\n"
        "    print(f'  {\"log10(B) [G]\":<22s} {med:.2f}  (+{hi-med:.2f}/-{med-lo:.2f})')\n"
        "    col = Pdot_ip_samples[mask]\n"
        "    med = np.median(col); lo, hi = np.percentile(col, [16, 84])\n"
        "    print(f'  {\"Pdot [s/s]\":<22s} {med:+.2e}  (+{hi-med:.2e}/{lo-med:.2e})')\n"
        "\n"
        "_branch_stats(acc_mask,  'Accretion')\n"
        "_branch_stats(prop_mask, 'Propeller')\n"
        "\n"
        "# --- Visualise regime split in the (M_WD, log B) plane ---\n"
        "fig, ax = plt.subplots(figsize=(6, 5))\n"
        "ax.scatter(samples_ip_eq[acc_mask, 0],  log10_B_samples_ip[acc_mask],\n"
        "           s=2, alpha=0.3, color='C0',\n"
        "           label=f'Accretion ({f_acc*100:.1f}%)')\n"
        "ax.scatter(samples_ip_eq[prop_mask, 0], log10_B_samples_ip[prop_mask],\n"
        "           s=2, alpha=0.3, color='C1',\n"
        "           label=f'Propeller ({f_prop*100:.1f}%)')\n"
        "ax.set_xlabel(r'$M_{\\rm WD}$ ($M_\\odot$)')\n"
        "ax.set_ylabel(r'$\\log_{10} B$ (G)')\n"
        "ax.set_title('IP posterior split by magnetospheric regime')\n"
        "ax.legend(markerscale=4, fontsize=10)\n"
        "fig.savefig(f'{savedir}/bayesian_ip_regime_split.pdf', bbox_inches='tight')\n"
        "plt.show()\n"
    ),
}

# Insert the new cells after cell index 12.
cells.insert(13, regime_md)
cells.insert(14, regime_code)

# IMPORTANT: after this insert the *original* cell indices shift by +2.
# Original -> new:
#   13 -> 15 (markdown "### Bayes Factor")
#   14 -> 16 (code:   Bayes factor computation)
#   15 -> 17 (markdown "### Corner Plots")
#   16 -> 18 (code:   IP corner plot)
#   17 -> 19 (code:   nmCV corner plot)
#   18 -> 20 (markdown "### Predicted Pdot")
#   19 -> 21 (code:   Pdot posterior hist)
#   20 -> 22 (markdown "### Posterior Predictive Check")
#   21 -> 23 (code:   PPC)
#   22 -> 24 (markdown "## Impact of Pdot Measurement")
#   23 -> 25 (code:   log_likelihood_ip_old + samplers)
#   24 -> 26 (code:   comparison table)
#   25 -> 27 (code:   IP overlay corner)
#   26 -> 28 (code:   nmCV overlay corner)
#   27 -> 29 (markdown summary)


# ---------------------------------------------------------------------------
# New cell index 18 — IP corner plot (originally cell 16): marginalise over xi
# ---------------------------------------------------------------------------
set_source(18, '''# --- IP corner plot (M_WD, log Mdot, log B) with prior overlays ---
# xi is marginalised by simply not displaying it (samples are already over xi).
import copy
from scipy.stats import gaussian_kde

results_ip_B = copy.deepcopy(results_ip)

# Rebuild the samples array with three display dimensions:
#   col 0: M_WD  (unchanged)
#   col 1: log10 Mdot  (unchanged)
#   col 2: log10 B  (derived from log10 mu at col 2 and R(M_WD))
# xi (original col 3) is dropped — marginalised.
raw = results_ip_B['samples']
R_arr = np.array([get_R(m) for m in raw[:, 0]])
logB_arr = np.log10(10**raw[:, 2] / R_arr**3)
display_samples = np.column_stack([raw[:, 0], raw[:, 1], logB_arr])
results_ip_B['samples'] = display_samples

# Compute B prior samples for overlay (log-uniform mu marginalised over N(M))
np.random.seed(42)
N_prior = 50000
M_prior = np.random.normal(POP_M_MU, POP_M_SIG, N_prior)
M_prior = np.clip(M_prior, M_min_table + 0.001, M_max_table - 0.001)
logmu_prior = np.random.uniform(LOGMU_LO, LOGMU_HI, N_prior)
R_prior = np.array([get_R(m) for m in M_prior])
logB_prior = logmu_prior - 3 * np.log10(R_prior)


def overlay_prior(ax, grid, pdf, **kwargs):
    """Scale prior PDF to match histogram area for correct normalization."""
    hist_area = sum(p.get_width() * p.get_height()
                    for p in ax.patches if p.get_height() > 0)
    if hist_area <= 0:
        return
    prior_area = np.trapezoid(pdf, grid)
    if prior_area <= 0:
        return
    pdf_scaled = pdf * (hist_area / prior_area)
    ax.plot(grid, pdf_scaled, **kwargs)
    ymax = max(ax.get_ylim()[1], pdf_scaled.max() * 1.05)
    ax.set_ylim([0, ymax])


# --- Prior-range view ---
span_ip_wide = [
    (POP_M_MU - 4*POP_M_SIG, POP_M_MU + 4*POP_M_SIG),
    (LOGMDOT_LO - 0.5, LOGMDOT_HI + 0.5),
    (np.min(logB_prior), np.max(logB_prior)-3),
]

plt.rcParams['figure.constrained_layout.use'] = False
fig, axes = dyplot.cornerplot(
    results_ip_B, span=span_ip_wide,
    labels=[r'$M_{\\rm WD}$ ($M_\\odot$)', r'$\\log_{10}\\dot{M}$ (g/s)', r'$\\log_{10} B$ (G)'],
    quantiles_2d=[0.1, 0.55, 0.93, 0.999],
    quantiles=None, show_titles=True,
    title_quantiles=[0.16, 0.5, 0.84], title_fmt='.3f',
    title_kwargs={'fontsize': 12}, label_kwargs={'fontsize': 13},
    color='C0',
)

# Overlay priors on diagonals
ax = axes[0, 0]
M_grid = np.linspace(*span_ip_wide[0], 300)
overlay_prior(ax, M_grid, norm.pdf(M_grid, POP_M_MU, POP_M_SIG),
              color='k', ls='--', lw=1.5, label='Prior')
ax = axes[1, 1]
Mdot_grid = np.linspace(*span_ip_wide[1], 300)
Mdot_prior_pdf = np.where((Mdot_grid >= LOGMDOT_LO) & (Mdot_grid <= LOGMDOT_HI),
                           1.0 / (LOGMDOT_HI - LOGMDOT_LO), 0.0)
overlay_prior(ax, Mdot_grid, Mdot_prior_pdf, color='k', ls='--', lw=1.5)
ax = axes[2, 2]
kde_B = gaussian_kde(logB_prior)
B_grid = np.linspace(*span_ip_wide[2], 300)
overlay_prior(ax, B_grid, kde_B(B_grid), color='k', ls='--', lw=1.5)

fig.suptitle(r'IP Model — prior range (marginalised over $\\xi$, with $L_X$ + $kT$)', fontsize=14, y=1.02)
fig.savefig(f'{savedir}/bayesian_corner_ip_spectral.pdf', bbox_inches='tight')
plt.show()
plt.rcParams['figure.constrained_layout.use'] = True

# --- Posterior zoom ---
N_SIGMA_ZOOM = 6

span_ip_zoom = []
for col in [samples_ip_eq[:, 0], samples_ip_eq[:, 1], log10_B_samples_ip]:
    med = np.median(col)
    sig = (np.percentile(col, 84) - np.percentile(col, 16)) / 2
    span_ip_zoom.append((med - N_SIGMA_ZOOM*sig, med + N_SIGMA_ZOOM*sig))

plt.rcParams['figure.constrained_layout.use'] = False
fig2, axes2 = dyplot.cornerplot(
    results_ip_B, span=span_ip_zoom,
    labels=[r'$M_{\\rm WD}$ ($M_\\odot$)', r'$\\log_{10}\\dot{M}$ (g/s)', r'$\\log_{10} B$ (G)'],
    quantiles=None, show_titles=True,
    title_quantiles=[0.16, 0.5, 0.84], title_fmt='.3f',
    title_kwargs={'fontsize': 12}, label_kwargs={'fontsize': 13},
    color='C0',
)
fig2.suptitle(r'IP Model — posterior zoom ($\\pm$' + f'{N_SIGMA_ZOOM}' + r'$\\sigma$, with $L_X$ + $kT$)',
              fontsize=14, y=1.02)
fig2.savefig(f'{savedir}/bayesian_corner_ip_spectral_zoom.pdf', bbox_inches='tight')
plt.show()
plt.rcParams['figure.constrained_layout.use'] = True
''')


# ---------------------------------------------------------------------------
# New cell index 21 — predicted Pdot distribution (originally cell 19).
# Pass xi into pdot_ip; reuse Pdot_ip_samples computed earlier if present.
# ---------------------------------------------------------------------------
set_source(21, '''# Compute Pdot from posterior samples
Pdot_nmcv_samples = np.array([
    pdot_nmcv(s[0], s[1]) for s in samples_nmcv_eq
])
# IP Pdot samples are already computed in the regime-partition cell as Pdot_ip_samples.
# Re-compute here only if that variable isn't available (e.g. cells run out of order).
try:
    Pdot_ip_samples
except NameError:
    Pdot_ip_samples = np.array([
        pdot_ip(s[0], s[1], s[2], s[3]) for s in samples_ip_eq
    ])

# --- Print summary statistics ---
for name, samp in [('nmCV', Pdot_nmcv_samples), ('IP', Pdot_ip_samples)]:
    med = np.median(samp)
    lo, hi = np.percentile(samp, [16, 84])
    mn, mx = samp.min(), samp.max()
    print(f"{name}:  Pdot = {med:.2e}  (+{hi-med:.2e} / {lo-med:.2e})")
    print(f"       range = [{mn:.2e}, {mx:.2e}]")
    frac_consistent = np.mean(np.abs(samp) < sigma_Pdot) * 100
    frac_2sig = np.mean(np.abs(samp) < 2*sigma_Pdot) * 100
    print(f"       {frac_consistent:.1f}% within 1-sigma of Pdot=0;  {frac_2sig:.1f}% within 2-sigma")
    print()

fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(7, 5), sharex=True,
                                      gridspec_kw={'height_ratios': [1, 1]})

bins = np.linspace(-5e-15, 5e-15, 80)
ax_top.hist(Pdot_nmcv_samples, bins=bins, density=True, alpha=0.6, color='C3', label='nmCV')
ax_bot.hist(Pdot_ip_samples, bins=bins, density=True, alpha=0.6, color='C0', label='IP')

for ax in [ax_top, ax_bot]:
    ax.axvline(0, color='gray', ls=':', lw=1)
    for nsig, ls, alph in [(1, '-', 0.5), (2, '--', 0.3)]:
        ax.axvline(nsig * sigma_Pdot, color='k', ls=ls, alpha=alph, lw=1.2)
        ax.axvline(-nsig * sigma_Pdot, color='k', ls=ls, alpha=alph, lw=1.2)
    ax.legend(fontsize=11)
    ax.set_ylabel('Density')

ax_bot.set_xlabel(r'$\\dot{P}$ (s/s)')
fig.suptitle(r'Predicted $\\dot{P}$ distributions (with $L_X$ + $kT$)', fontsize=13)
fig.savefig(f'{savedir}/bayesian_pdot_posterior_spectral.pdf', bbox_inches='tight')
plt.show()
''')


# ---------------------------------------------------------------------------
# New cell index 23 — posterior predictive check (originally cell 21).
# Pass xi into kT_ip / LX_ip.
# ---------------------------------------------------------------------------
set_source(23, '''# --- Predicted kT and L_X from posterior ---
kT_nmcv_pred = np.array([kT_nmcv(s[0]) for s in samples_nmcv_eq])
LX_nmcv_pred = np.array([LX_nmcv(s[0], s[1]) for s in samples_nmcv_eq])

kT_ip_pred = np.array([kT_ip(s[0], s[1], s[2], s[3]) for s in samples_ip_eq])
LX_ip_pred = np.array([LX_ip(s[0], s[1], s[2], s[3]) for s in samples_ip_eq])
# Accretion-branch samples only have defined spectral predictions.
# Propeller-branch samples yield NaN from kT_ip/LX_ip — filter them out here.
valid_ip = np.isfinite(kT_ip_pred) & np.isfinite(LX_ip_pred)

# --- Chi-squared breakdown (accretion-branch IP only for kT/L_X) ---
print("=== Per-observable chi^2 (median of posterior) ===\\n")
for name, Pdot_s, kT_s, LX_s in [
    ('nmCV', Pdot_nmcv_samples, kT_nmcv_pred, LX_nmcv_pred),
    ('IP (accretion branch)', Pdot_ip_samples[valid_ip],
     kT_ip_pred[valid_ip], LX_ip_pred[valid_ip]),
]:
    chi2_Pdot = np.median(((Pdot_obs - Pdot_s) / sigma_Pdot)**2)
    chi2_kT = np.median(((kT_obs - kT_s) / sigma_kT)**2)
    chi2_LX = np.median(((L_X_obs - LX_s) / sigma_LX)**2)
    print(f"{name}:")
    print(f"  chi2(Pdot) = {chi2_Pdot:.2f}")
    print(f"  chi2(kT)   = {chi2_kT:.2f}")
    print(f"  chi2(L_X)  = {chi2_LX:.2f}")
    print(f"  total      = {chi2_Pdot + chi2_kT + chi2_LX:.2f}")
    print()

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# kT panel
ax = axes[0]
ax.hist(kT_nmcv_pred, bins=50, density=True, alpha=0.5, color='C3', label='nmCV')
ax.hist(kT_ip_pred[valid_ip], bins=50, density=True, alpha=0.5, color='C0',
        label='IP (accretion branch)')
ax.axvline(kT_obs, color='k', lw=2, label=f'Observed ({kT_obs} keV)')
ax.axvspan(kT_obs - sigma_kT, kT_obs + sigma_kT, alpha=0.15, color='gray')
ax.set_xlabel(r'$kT_{\\rm max}$ (keV)')
ax.set_ylabel('Density')
ax.legend(fontsize=10)

# L_X panel
ax = axes[1]
ax.hist(LX_nmcv_pred, bins=50, density=True, alpha=0.5, color='C3', label='nmCV')
ax.hist(LX_ip_pred[valid_ip], bins=50, density=True, alpha=0.5, color='C0',
        label='IP (accretion branch)')
ax.axvline(L_X_obs, color='k', lw=2, label=f'Observed')
ax.axvspan(L_X_obs - sigma_LX, L_X_obs + sigma_LX, alpha=0.15, color='gray')
ax.set_xlabel(r'$L_X$ (erg/s)')
ax.set_ylabel('Density')
ax.legend(fontsize=10)

fig.suptitle('Posterior Predictive Check', fontsize=13)
fig.savefig(f'{savedir}/bayesian_predictive_check_spectral.pdf', bbox_inches='tight')
plt.show()
''')


# ---------------------------------------------------------------------------
# New cell index 25 — old-Pdot re-run (originally cell 23). Update to 4D IP.
# ---------------------------------------------------------------------------
set_source(25, '''# --- Re-run with old (weak) Pdot constraint ---
# Same sampler + priors + branch-aware likelihood as the main runs so that
# the ln Z comparison between "new" and "old" Pdot constraints is
# apples-to-apples for the Bayes factor.
sigma_Pdot_old = 1.8e-12  # old 1-sigma: essentially no constraint

def log_likelihood_nmcv_old(theta):
    M_wd, log_Mdot = theta
    if M_wd <= M_min_table or M_wd >= M_max_table:
        return -1e100
    Pdot_pred = pdot_nmcv(M_wd, log_Mdot)
    kT_pred = kT_nmcv(M_wd)
    LX_pred = LX_nmcv(M_wd, log_Mdot)
    chi2  = ((Pdot_obs - Pdot_pred) / sigma_Pdot_old)**2
    chi2 += ((kT_obs - kT_pred) / sigma_kT)**2
    chi2 += ((L_X_obs - LX_pred) / sigma_LX)**2
    return -0.5 * chi2


def log_likelihood_ip_old(theta):
    M_wd, log_Mdot, log_mu, xi_s = theta
    if M_wd <= M_min_table or M_wd >= M_max_table:
        return -1e100

    Pdot_pred = pdot_ip(M_wd, log_Mdot, log_mu, xi_s)
    if np.isnan(Pdot_pred):
        return -1e100

    chi2 = ((Pdot_obs - Pdot_pred) / sigma_Pdot_old)**2

    M_cgs = M_wd * Msun
    Mdot = 10**log_Mdot
    mu = 10**log_mu
    r_m = xi_s * (mu**4 / (2 * G * M_cgs * Mdot**2))**(1.0/7.0)
    r_co = (G * M_cgs * P**2 / (4 * np.pi**2))**(1.0/3.0)

    if r_m <= r_co:
        kT_pred = kT_ip(M_wd, log_Mdot, log_mu, xi_s)
        LX_pred = LX_ip(M_wd, log_Mdot, log_mu, xi_s)
        if np.isnan(kT_pred) or kT_pred <= 0:
            return -1e100
        if np.isnan(LX_pred) or LX_pred <= 0:
            return -1e100
        chi2 += ((kT_obs - kT_pred) / sigma_kT)**2
        chi2 += ((L_X_obs - LX_pred) / sigma_LX)**2
    else:
        LX_resid = eta * G * M_cgs * Mdot / r_m
        LX_cap = L_X_obs + 2 * sigma_LX
        if LX_resid > LX_cap:
            chi2 += ((LX_resid - LX_cap) / sigma_LX)**2

    return -0.5 * chi2


print("Running IP sampling with OLD Pdot constraint (4D)...")
sampler_ip_old = dynesty.DynamicNestedSampler(
    log_likelihood_ip_old, prior_transform_ip, ndim=4,
    nlive=2000, bound='multi', sample='rwalk'
)
sampler_ip_old.run_nested(print_progress=True)
results_ip_old = sampler_ip_old.results
print(f"IP (old): logZ = {results_ip_old['logz'][-1]:.2f} +/- {results_ip_old['logzerr'][-1]:.2f}")

print("\\nRunning nmCV sampling with OLD Pdot constraint...")
sampler_nmcv_old = dynesty.DynamicNestedSampler(
    log_likelihood_nmcv_old, prior_transform_nmcv, ndim=2,
    nlive=2000, bound='multi', sample='rwalk'
)
sampler_nmcv_old.run_nested(print_progress=True)
results_nmcv_old = sampler_nmcv_old.results
print(f"nmCV (old): logZ = {results_nmcv_old['logz'][-1]:.2f} +/- {results_nmcv_old['logzerr'][-1]:.2f}")
''')


# ---------------------------------------------------------------------------
# New cell index 26 — comparison table (originally cell 24).
# Add xi row and keep existing M/Mdot/mu/B reporting.
# ---------------------------------------------------------------------------
set_source(26, '''# --- Extract old-constraint posteriors ---
w_ip_old = np.exp(results_ip_old['logwt'] - results_ip_old['logz'][-1])
w_ip_old /= w_ip_old.sum()
samples_ip_old_eq = dyfunc.resample_equal(results_ip_old['samples'], w_ip_old)

w_nmcv_old = np.exp(results_nmcv_old['logwt'] - results_nmcv_old['logz'][-1])
w_nmcv_old /= w_nmcv_old.sum()
samples_nmcv_old_eq = dyfunc.resample_equal(results_nmcv_old['samples'], w_nmcv_old)

# Derived B for IP old
R_old = np.array([get_R(m) for m in samples_ip_old_eq[:, 0]])
log10_B_old = samples_ip_old_eq[:, 2] - 3 * np.log10(R_old)

# --- Comparison table ---
def fmt(samples, idx):
    med = np.median(samples[:, idx])
    lo, hi = np.percentile(samples[:, idx], [16, 84])
    return f"{med:.3f} (+{hi-med:.3f}/-{med-lo:.3f})"

print("=== Posterior comparison: new vs old Pdot constraint ===\\n")
header = f"{'Parameter':<22s} {'nmCV (new)':<28s} {'nmCV (old)':<28s}"
print(header)
print("-" * len(header))
print(f"{'M_WD [Msun]':<22s} {fmt(samples_nmcv_eq, 0):<28s} {fmt(samples_nmcv_old_eq, 0):<28s}")
print(f"{'log10(Mdot) [g/s]':<22s} {fmt(samples_nmcv_eq, 1):<28s} {fmt(samples_nmcv_old_eq, 1):<28s}")
print()

header = f"{'Parameter':<22s} {'IP (new)':<28s} {'IP (old)':<28s}"
print(header)
print("-" * len(header))
print(f"{'M_WD [Msun]':<22s} {fmt(samples_ip_eq, 0):<28s} {fmt(samples_ip_old_eq, 0):<28s}")
print(f"{'log10(Mdot) [g/s]':<22s} {fmt(samples_ip_eq, 1):<28s} {fmt(samples_ip_old_eq, 1):<28s}")
print(f"{'log10(mu) [G cm^3]':<22s} {fmt(samples_ip_eq, 2):<28s} {fmt(samples_ip_old_eq, 2):<28s}")
print(f"{'xi':<22s} {fmt(samples_ip_eq, 3):<28s} {fmt(samples_ip_old_eq, 3):<28s}")

med_B_new = np.median(log10_B_samples_ip)
med_B_old = np.median(log10_B_old)
lo_B_new, hi_B_new = np.percentile(log10_B_samples_ip, [16, 84])
lo_B_old, hi_B_old = np.percentile(log10_B_old, [16, 84])
new_B = f"{med_B_new:.2f} (+{hi_B_new-med_B_new:.2f}/-{med_B_new-lo_B_new:.2f})"
old_B = f"{med_B_old:.2f} (+{hi_B_old-med_B_old:.2f}/-{med_B_old-lo_B_old:.2f})"
print(f"{'log10(B) [G]':<22s} {new_B:<28s} {old_B:<28s}")

# --- Regime fraction comparison ---
def _regime_frac(samples_eq):
    M_cgs = samples_eq[:, 0] * Msun
    Mdot = 10**samples_eq[:, 1]
    mu_ = 10**samples_eq[:, 2]
    xi_ = samples_eq[:, 3]
    r_m = xi_ * (mu_**4 / (2 * G * M_cgs * Mdot**2))**(1.0/7.0)
    r_co = (G * M_cgs * P**2 / (4 * np.pi**2))**(1.0/3.0)
    return (r_m <= r_co).mean(), (r_m > r_co).mean()

f_acc_new, f_prop_new = _regime_frac(samples_ip_eq)
f_acc_old, f_prop_old = _regime_frac(samples_ip_old_eq)
print(f"\\n{'Regime':<22s} {'IP (new)':<28s} {'IP (old)':<28s}")
print("-" * len(header))
print(f"{'accretion fraction':<22s} {f_acc_new:<28.3f} {f_acc_old:<28.3f}")
print(f"{'propeller fraction':<22s} {f_prop_new:<28.3f} {f_prop_old:<28.3f}")

print("\\n=== Bayes factor comparison ===")
logZ_ip_new = results_ip['logz'][-1]
logZ_nmcv_new = results_nmcv['logz'][-1]
logZ_ip_o = results_ip_old['logz'][-1]
logZ_nmcv_o = results_nmcv_old['logz'][-1]
lnBF_new = logZ_ip_new - logZ_nmcv_new
lnBF_old = logZ_ip_o - logZ_nmcv_o
print(f"  ln(BF_IP/nmCV) new Pdot: {lnBF_new:+.2f}  (BF = {np.exp(lnBF_new):.4f})")
print(f"  ln(BF_IP/nmCV) old Pdot: {lnBF_old:+.2f}  (BF = {np.exp(lnBF_old):.4f})")
''')


# ---------------------------------------------------------------------------
# New cell index 27 — IP overlay corner (originally cell 25). Marginalise xi.
# ---------------------------------------------------------------------------
set_source(27, '''# --- Overlay IP corner plots: new (blue) vs old (orange) Pdot constraint ---
# Both results are 4D; drop the xi column and replace log10(mu) with log10(B).
results_ip_B_new = copy.deepcopy(results_ip)
raw_new = results_ip_B_new['samples']
R_tmp_new = np.array([get_R(m) for m in raw_new[:, 0]])
logB_new = np.log10(10**raw_new[:, 2] / R_tmp_new**3)
results_ip_B_new['samples'] = np.column_stack([raw_new[:, 0], raw_new[:, 1], logB_new])

results_ip_B_old = copy.deepcopy(results_ip_old)
raw_old = results_ip_B_old['samples']
R_tmp_old = np.array([get_R(m) for m in raw_old[:, 0]])
logB_old_col = np.log10(10**raw_old[:, 2] / R_tmp_old**3)
results_ip_B_old['samples'] = np.column_stack([raw_old[:, 0], raw_old[:, 1], logB_old_col])

# Span encompassing both posteriors
all_ip = np.vstack([samples_ip_eq, samples_ip_old_eq])
all_B = np.concatenate([log10_B_samples_ip, log10_B_old])
span_ip_cmp = []
for col in [all_ip[:, 0], all_ip[:, 1], all_B]:
    lo, hi = np.percentile(col, [0.5, 99.5])
    w = hi - lo
    span_ip_cmp.append((lo - 0.3*w, hi + 0.3*w))

plt.rcParams['figure.constrained_layout.use'] = False
fig, axes = dyplot.cornerplot(
    results_ip_B_new, span=span_ip_cmp,
    labels=[r'$M_{\\rm WD}$ ($M_\\odot$)', r'$\\log_{10}\\dot{M}$ (g/s)', r'$\\log_{10} B$ (G)'],
    quantiles=None, show_titles=True, title_quantiles=[0.16, 0.5, 0.84],
    title_fmt='.3f', title_kwargs={'fontsize': 12}, label_kwargs={'fontsize': 13},
    color='C0',
)
dyplot.cornerplot(
    results_ip_B_old, span=span_ip_cmp,
    labels=[r'$M_{\\rm WD}$ ($M_\\odot$)', r'$\\log_{10}\\dot{M}$ (g/s)', r'$\\log_{10} B$ (G)'],
    quantiles=None, show_titles=False,
    title_kwargs={'fontsize': 12}, label_kwargs={'fontsize': 13},
    color='C1', fig=(fig, axes),
)
from matplotlib.lines import Line2D
fig.legend(handles=[
    Line2D([0], [0], color='C0', lw=2, label=r'New $\\sigma_{\\dot{P}}$ = 1.7e-15'),
    Line2D([0], [0], color='C1', lw=2, label=r'Old $\\sigma_{\\dot{P}}$ = 1.8e-12'),
], loc='upper right', fontsize=11, frameon=True)
fig.suptitle(r'IP Model: $\\dot P$ constraint comparison (marginalised over $\\xi$)', fontsize=14, y=1.02)
fig.savefig(f'{savedir}/bayesian_corner_ip_pdot_comparison.pdf', bbox_inches='tight')
plt.show()
plt.rcParams['figure.constrained_layout.use'] = True
''')


# ---------------------------------------------------------------------------
# Write back
# ---------------------------------------------------------------------------
nb_path.write_text(json.dumps(nb, indent=1) + "\n")
print(f"Patched {nb_path.name} ({len(cells)} cells, backup at {backup.name})")
