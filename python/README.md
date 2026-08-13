# tilin
## PASCAL implementation of the Thiele–Innes technique for determining the orbits of binary stars 

This older PASCAL code was utilized in a series of papers published by astronomers from the Pulkovo Observatory. The author of this code is Igor Izmailov. A good instance of the paper based on the presented elaboration is  Izmailov, I.~S.\ 2019.\ The Orbits of 451 Wide Visual Double Stars.\ Astronomy Letters 45, 30–38. doi:10.1134/S106377371901002Xhttps://ui.adsabs.harvard.edu/abs/2019AstL...45...30I/abstract. The code is based on the Thiele–Innes method and allows for calculating the Campbell elements of the relative orbit for the binary star. 

### Compilation 
The code can be compiled using Turbo Pascal on any Windows OS (see tilin132.exe in the file set) or with the FreePascal compiler (Linux, macOS, etc.). The root file is TILIN132.PAS. The command in the Linux terminal appears as follows:
```
fpc -Mtp TILIN132.PAS
```
Of course, you should be in the root directory of the project. The build procedure result is the TILIN132 file. 

### Utilization

It is easy to execute:
```
./TILIN132
```
The input.txt file contains a table with relative positions, including separations in arcsec (\rho) and positional angles in degrees (\theta) for the corresponding epochs (in years).

Here is an example of the input.txt:
```
epoch (yr) \theta (deg) \rho (arcsec) 
1836.21000 295.60397   2.50000
1852.92000 292.50357   2.89000
1857.90000 292.10344   2.63000
1880.59000 295.70289   2.10000
1891.84000 293.00262   2.04000
```
You can use your input.txt with your data according to the presented form.

The calculation results, including the Campbell elements, will be saved in output.txt (see an example below).
```
t0            =  1949.885285185185185
rho0          =     2.053027044710415
theta0        =   310.995348392656718
rho t         =    -0.007180785989696
theta t       =     0.254708783231639
rad.curvature =     1.844834589269195
R V["/year]   =    -0.078449413947353
mass(pi=1)    =     0.000028641534337

a    =     3.096683902242133
e    =     0.832648372013274
T    =  2058.850225450467577
P    =  1018.233016244432854
omega=   263.998907095081453
Node =   175.906109257159364
i    =    55.125493329908873

RMS=     0.139600983918518

RMS
all obs.=     0.139600983918518

    Ep       Rho(obs)      Theta(obs)     Rho(c)     Theta(c)       dRho        dTheta          X            Y 
1836.21000    2.5000000  295.6039700    2.6507819  289.1277343   -0.1507819    6.4762357   -2.2545065    1.0803706
1852.92000    2.8900000  292.5035700    2.5823538  291.7476326    0.3076462    0.7559374   -2.6699429    1.1061215
1857.90000    2.6300000  292.1034400    2.5609847  292.5560334    0.0690153   -0.4525934   -2.4367109    0.9896161
1880.59000    2.1000000  295.7028900    2.4573199  296.4240501   -0.3573199   -0.7211601   -1.8922158    0.9107795
```
Rho(obs), Theta(obs) are taken from input.txt.
Rho(c), Theta(c) are calculated using estimated a, e, T, P, omega, Node, and i.
dRho,dTheta are O minus C (i.e. Rho(obs) - Rho(c) and Theta(obs)-Theta(c)).

### Python module

`tilin_module.py` wraps the same fitting code for use from other Python
scripts, in place of `input.txt`/`output.txt`. Observations are passed in as
a pandas DataFrame (with configurable column names for epoch, separation,
and position angle), and results come back as a dataclass instead of a file:

```python
import pandas as pd
from tilin_module import fit_orbit_dataframe

obs = pd.DataFrame({
    'year':  [1836.21, 1852.92, 1857.90, 1880.59, 1891.84],
    'rho':   [2.50,    2.89,    2.63,    2.10,    2.04],
    'theta': [295.60,  292.50,  292.10,  295.70,  293.00],
})

result = fit_orbit_dataframe(obs, year_col='year', rho_col='rho', theta_col='theta')

print(result.elements)    # {'a':..., 'e':..., 'T':..., 'P':..., 'omega':..., 'Node':..., 'i':...}
print(result.rms)         # weighted RMS residual (arcsec) over the used observations
print(result.residuals)   # per-observation O-C table, indexed like `obs`, with an 'excluded' flag
```

Pass `verbose=False` to suppress the fit-progress messages printed to stdout.
Rows with NaN in the year/rho/theta columns are dropped before fitting.

### Orbital-element uncertainty (bootstrap)

`tilin_module.py` also has a parametric-bootstrap estimator for the fitted
elements' uncertainty. It takes the fitted model's predicted sky position at
each observation epoch, adds isotropic Gaussian noise (std = the fit's RMS,
in arcsec) to build synthetic observation sets, refits each one from scratch,
and collects the resulting elements:

```python
from tilin_module import bootstrap_orbit_elements, orbit_element_limits, plot_orbit_corner

samples = bootstrap_orbit_elements(result, n_boot=100)   # ~100x the cost of one fit
limits  = orbit_element_limits(samples)                  # 16th/84th-percentile (~1σ) bounds per element
plot_orbit_corner(samples, best_fit=result.elements, filename='orbit_corner.png')
```

`samples` is a DataFrame with one row per successful bootstrap fit
(`a, e, T, P, omega, Node, i, rms`). `orbit_element_limits` returns a
`lower`/`upper` bound per element (percentiles are configurable via
`lower_pct`/`upper_pct`). `plot_orbit_corner` saves a PNG corner plot
(pairwise scatter + histograms, with the nominal fit marked in red) — useful
for spotting degenerate/multimodal solutions on short observation arcs,
which is common for wide visual binaries. Note that `omega`/`Node` are
angles that can wrap around 0/360°; a distribution straddling the wrap point
will make the plain percentile bounds misleading, so check the plot.

Since each bootstrap realization reruns the full fit (global search +
outlier rejection), `n_boot=100` takes roughly 100x as long as a single
`fit_orbit_dataframe` call — seconds to a few minutes depending on the
number of observations. `bootstrap_orbit_elements` shows a live progress
bar by default (via `tqdm.auto` — a plain bar in a terminal, a widget in
Jupyter; pass `progress=False` to turn it off). This needs the optional
`tqdm` package; without it the bar is silently skipped.

`plot_orbit_family(samples, result, filename='orbit_family.png')` saves a
complementary view: the bootstrap orbits drawn as curves in the plane of the
sky (X/Y, arcsec), the nominal fit in red, and the observations as points
(hollow for any excluded as outliers). This is often more intuitive than the
corner plot for spotting *how* the solutions diverge — e.g. several
differently-sized loops that all thread the same short observed arc. The
default view is zoomed to the nominal orbit + observations (pass
`zoom_pad=None` to instead autoscale to the full bootstrap family, including
any wide/degenerate outlier orbits).
