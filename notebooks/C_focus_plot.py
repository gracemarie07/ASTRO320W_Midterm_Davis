"""Draw the focus measurements listed in the report's AstroImageJ table."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
focus = np.array([17000,17000,17000,16500,16000,16000,15500,15000])
fwhm = np.array([11.55,12.41,11.19,10.30,9.75,10.68,10.72,12.74])
x = np.linspace(15000,17000,200)
fit = np.polyfit(focus-16000,fwhm,2)
plt.rcParams.update({'font.size':11})
fig, ax = plt.subplots(figsize=(6.5,3.6),layout='constrained')
ax.plot(x,np.polyval(fit,x-16000),color='black',lw=1)
ax.plot(focus,fwhm,'o',color='#245884',ms=4)
for position in np.unique(focus):
    values = fwhm[focus==position]
    if len(values)>1:
        ax.errorbar(position,values.mean(),yerr=values.std(ddof=1),fmt='s',color='black',capsize=3,ms=3)
ax.set(xlabel='Focus position (micrometers)',ylabel='FWHM (pixels)')
ax.grid(alpha=.2)
for folder in ['report','results']:
    fig.savefig(ROOT/folder/'focus_panels.pdf',bbox_inches='tight')
plt.close(fig)
