"""Draw the Ogma exposure-time figure from the original AstroImageJ export."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/astro-excel-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
measurements = pd.read_csv(ROOT/'data/newtable.csv')
measurements = measurements.loc[measurements.Source_Radius == 30].sort_values('EXPTIME')
assert len(measurements) == 12
time = measurements.EXPTIME.to_numpy()
counts = measurements['Source-Sky_T1'].to_numpy()
width = measurements.FWHM_T1.to_numpy()
plt.rcParams.update({'font.size':11, 'axes.labelsize':11, 'axes.titlesize':11,
                     'xtick.labelsize':11, 'ytick.labelsize':11})
fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.8), layout='constrained')
for ax, values, title, ylabel in zip(axes, [width, counts/1e6, counts/time/1e5],
        ['(a) Image width', '(b) Net counts', '(c) Count rate'],
        ['FWHM (pixels)', 'Net counts (million DN)', 'Rate (100000 DN/s)']):
    ax.plot(time, values, 'o-', color='#245884', markersize=3, linewidth=.8)
    ax.set(title=title, xlabel='Time (s)', ylabel=ylabel, xticks=[0,150,300])
    ax.grid(alpha=.2)
fig.savefig(ROOT/'report/exposure_panels.pdf', bbox_inches='tight')
fig.savefig(ROOT/'results/exposure_panels.pdf', bbox_inches='tight')
fig.savefig(ROOT/'results/exposure_panels.png', dpi=220, bbox_inches='tight')
plt.close(fig)
print('Ogma width, net counts and count-rate panels saved from the twelve original measurements.')
