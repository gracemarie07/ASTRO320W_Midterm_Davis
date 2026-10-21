"""Reproduce representative flat comparisons from explicitly identified raw frames."""
from pathlib import Path
import os,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[1]
BASE=Path(os.environ.get('ASTRO320W_ORGANIZED',str(ROOT/'organized_data')))
files=[BASE/'B/Sky_Flats_R/Autosave Image -0002skyflat.fit',BASE/'B/Dome_Paired/DOME_FLAT_20_EXP-000110.fit']
for p in files:
 if not p.exists():raise FileNotFoundError(f'Set ASTRO320W_ORGANIZED; missing {p}')
plt.rcParams.update({'font.size':11.5,'axes.labelsize':11.5,'axes.titlesize':11.5,'xtick.labelsize':11.5,'ytick.labelsize':11.5})
f,axs=plt.subplots(2,2,figsize=(6.5,4.8),layout='constrained');stats=[]
for i,(p,name) in enumerate(zip(files,['Sky','Dome'])):
 a=fits.getdata(p).astype(float);h=fits.getheader(p)
 # Dome flats were dark-subtracted on acquisition; remove the signed stored pedestal once.
 if 'D' in str(h.get('CALSTAT','')):a+=float(h.get('PEDESTAL',0))
 y,x=np.array(a.shape)//2;hp=np.median(a[y-5:y+6,:],axis=0);vp=np.median(a[:,x-5:x+6],axis=1);median=np.median(a)
 for j,(v,label) in enumerate(zip([hp,vp],['Horizontal','Vertical'])):
  axs[i,j].plot(np.arange(len(v)),v/median,lw=.6);axs[i,j].set(xlabel='Pixel position',ylabel='Signal / frame median',title=f'{name}: {label.lower()}');axs[i,j].grid(alpha=.2)
 stats.append({'file':p.name,'exposure_s':h.get('EXPTIME'),'horizontal_min_max_normalized':[float(hp.min()/median),float(hp.max()/median)],'vertical_min_max_normalized':[float(vp.min()/median),float(vp.max()/median)]})
f.savefig(ROOT/'report/flat_profile_panels.pdf',bbox_inches='tight');plt.close(f)
shutil.copy2(ROOT/'report/flat_profile_panels.pdf',ROOT/'results/flat_profile_panels.pdf')
(ROOT/'results/flat_profile_provenance.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats,indent=2))
f,axs=plt.subplots(1,2,figsize=(6.5,3.7),layout='constrained')
for ax,p,title in zip(axs,files,['Sky: 6 s','Dome: 10 s']):
 a=fits.getdata(p).astype(float);h=fits.getheader(p)
 if 'D' in str(h.get('CALSTAT','')):a+=float(h.get('PEDESTAL',0))
 lo,hi=np.percentile(a,[3,99]);z=ax.imshow(a,origin='upper',cmap='gray',vmin=lo,vmax=hi);ax.set(title=title,xlabel='Array x (pixels)',ylabel='Array y (pixels)');f.colorbar(z,ax=ax,orientation='horizontal',label='DN',shrink=.9)
f.savefig(ROOT/'report/flat_examples_identified.pdf',bbox_inches='tight');plt.close(f)
shutil.copy2(ROOT/'report/flat_examples_identified.pdf',ROOT/'results/flat_examples_identified.pdf')
