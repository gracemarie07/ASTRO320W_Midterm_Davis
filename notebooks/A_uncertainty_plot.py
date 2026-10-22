"""Plot the Exercise A dark measurements and linear fits."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/astro-excel-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
time = pd.read_csv(ROOT / 'data/dark_time_with_errors.csv')
temp = pd.read_csv(ROOT / 'data/dark_temperature_with_errors.csv')

def linear_fit(x, y):
    # Same unweighted least-squares fit as the Excel workbook.
    coefficients, covariance = np.polyfit(x, y, 1, cov=True)
    return coefficients, covariance

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 11,
                     'axes.titlesize': 11, 'legend.fontsize': 10.5,
                     'xtick.labelsize': 11, 'ytick.labelsize': 11})
fig = plt.figure(figsize=(6.5, 7.4), layout='constrained')
grid_spec = fig.add_gridspec(2, 2)
axes = [fig.add_subplot(grid_spec[0, :]), fig.add_subplot(grid_spec[1, 0]),
        fig.add_subplot(grid_spec[1, 1])]
blue = '#245884'

ax = axes[0]
x, y = time.t_s.to_numpy(), time.Q_e.to_numpy()
grid = np.linspace(0, 60, 160)
time_fit, time_cov = linear_fit(x, y)
ax.plot(grid, np.polyval(time_fit, grid), color='black', lw=1, label='Linear fit')
ax.errorbar(x, y, yerr=time.Q_region_error_e,
            fmt='o', ms=2.5, color=blue, ecolor='black', capsize=4,
            elinewidth=1, label='Measurements')
ax.set(title='(a) Charge vs. exposure time', xlabel='Exposure time (s)',
       ylabel='Charge (electrons/pixel; log scale)')
ax.set_yscale('log')
ax.set_ylim(30, 1050)

ax = axes[1]
ax.errorbar(temp.T_start_C, temp.Q_e, yerr=temp.Q_region_error_e,
            fmt='o', ms=2.5, color=blue, ecolor='black', capsize=4)
ax.set(title='(b) Charge vs. temperature', xlabel='Start temperature (C)',
       ylabel='Charge (electrons/pixel; log scale)')

ax.set_yscale('log')
ax.set_ylim(5, 1500)

ax = axes[2]
inverse = 1000/(temp.T_start_C.to_numpy()+273.15)
log_rate = np.log(temp.I_e_s.to_numpy())
# Transform the rate's +/- one-standard-error endpoints. The cold point has a large fractional error.
rate, rate_error = temp.I_e_s.to_numpy(), temp.I_region_error_e_s.to_numpy()
assert np.all(rate > rate_error)
log_error = np.vstack([log_rate-np.log(rate-rate_error),
                       np.log(rate+rate_error)-log_rate])
grid = np.linspace(inverse.min(), inverse.max(), 160)
thermal_fit, thermal_cov = linear_fit(inverse, log_rate)
ax.plot(grid, np.polyval(thermal_fit, grid), color='black', lw=1)
ax.errorbar(inverse, log_rate, yerr=log_error, fmt='o', ms=2.5,
            color=blue, ecolor='black', capsize=4, elinewidth=1)
ax.set(title='(c) Linearized thermal fit', xlabel='1000 / T (1/K)',
       ylabel='ln(I / Iref)')

for ax in axes:
    ax.grid(alpha=.2)

fig.savefig(ROOT / 'report/dark_panels.pdf', bbox_inches='tight')
fig.savefig(ROOT / 'results/dark_panels.pdf', bbox_inches='tight')
fig.savefig(ROOT / 'results/dark_panels.png', dpi=220, bbox_inches='tight')
plt.close(fig)

print('Time-fit slope and intercept:', time_fit)
print('Their standard errors:', np.sqrt(np.diag(time_cov)))
