"""Measure bias maps and raw stellar peaks from FITS images."""
from pathlib import Path
import os, json, csv, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import linregress
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; OUT=ROOT/'results'; REPORT=ROOT/'report'
OUT.mkdir(exist_ok=True)
RAW=Path(os.environ.get('ASTRO320W_DATA',str(ROOT/'raw')))
BIAS=Path(os.environ.get('ASTRO320W_BIAS',str(RAW/'Exercise A'/'Bias_Frames')))
ASTRO=Path(os.environ.get('ASTRO320W_ORGANIZED',str(ROOT/'organized_data')))
plt.rcParams.update({'font.size':11,'axes.labelsize':11,'axes.titlesize':11,'legend.fontsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'savefig.dpi':200})
summary={}
def save(fig,name):
 fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');fig.savefig(OUT/(name+'.png'),bbox_inches='tight');plt.close(fig)
 shutil.copy2(OUT/(name+'.pdf'),REPORT/(name+'.pdf'))
def export(df,name):df.to_csv(DATA/name,index=False)
# Full-detector noise maps: disjoint-pair noise and mean single-frame spatial SD.
files=sorted(BIAS.glob('biasframe-*.fit'))
if len(files)==20:
 stack=np.stack([fits.getdata(f).astype(np.float32) for f in files]);master=stack.mean(axis=0);h,w=master.shape;grid=[];RN=np.zeros((8,8));SP=np.zeros((8,8))
 central=stack[:,1478:1578,1478:1578].mean(axis=(1,2),dtype=np.float64)
 export(pd.DataFrame({'file':[p.name for p in files],'region_pixels':10000,'regional_mean_DN':central}),'bias_region_means.csv')
 for iy,(y0,y1) in enumerate(zip(np.linspace(0,h,9,dtype=int)[:-1],np.linspace(0,h,9,dtype=int)[1:])):
  for ix,(x0,x1) in enumerate(zip(np.linspace(0,w,9,dtype=int)[:-1],np.linspace(0,w,9,dtype=int)[1:])):
   sub=stack[:,y0:y1,x0:x1];diff=sub[::2]-sub[1::2]
   RN[iy,ix]=np.mean(np.std(diff,axis=(1,2),ddof=1)/np.sqrt(2));SP[iy,ix]=np.mean(np.std(sub,axis=(1,2),ddof=1))
   grid.append([ix,iy,x0,x1,y0,y1,RN[iy,ix],SP[iy,ix]])
 export(pd.DataFrame(grid,columns=['block_x','block_y','x_start','x_stop','y_start','y_stop','paired_read_noise_DN','single_frame_spatial_SD_DN']),'bias_block_statistics.csv')
 fig,axs=plt.subplots(1,2,figsize=(6.5,3.5),layout='constrained')
 for ax,a,title in zip(axs,[RN,SP],['Paired read noise','Single-frame bias scatter']):
  z=ax.imshow(a,origin='upper',cmap='viridis');fig.colorbar(z,ax=ax,label='DN',shrink=.75);ax.set(title=title,xlabel='Block x',ylabel='Block y')
 save(fig,'bias_block_maps');summary['bias_maps']={'read_noise_min':float(RN.min()),'read_noise_max':float(RN.max()),'spatial_min':float(SP.min()),'spatial_max':float(SP.max())}

 del stack
# Raw stellar peak/plateau checks, without equating a sub-65535 plateau to ADC clipping.
focus=[]
for i in [1,2,3,4,6,7,8,9]:
 f=RAW/'EXERCISE C-2'/f'CCD Image {i}.fit'
 if not f.exists():continue
 im=fits.getdata(f).astype(float);y,x=np.unravel_index(np.argmax(im),im.shape);box=im[max(y-20,0):y+21,max(x-20,0):x+21];peak=box.max();focus.append(dict(image=i,peak_DN=peak,pixels_above_48000=int((box>48000).sum()),x=x,y=y))
 if i==6:
  fig,axs=plt.subplots(1,2,figsize=(6.5,3.0),layout='constrained');axs[0].imshow(box,origin='upper',cmap='gray');axs[0].set(title='Raw 6-s focus image',xlabel='Cutout x (pixels)',ylabel='Cutout y (pixels)');axs[1].plot(np.arange(-20,21),box[20,:]);axs[1].set(xlabel='Offset (pixels)',ylabel='Raw signal (DN)',title='Row through peak');save(fig,'focus_raw_check')
if focus:
 export(pd.DataFrame(focus),'focus_saturation_check.csv');summary['focus_check']={'min_peak':min(r['peak_DN'] for r in focus),'max_peak':max(r['peak_DN'] for r in focus),'min_plateau_pixels':min(r['pixels_above_48000'] for r in focus),'max_plateau_pixels':max(r['pixels_above_48000'] for r in focus)}
long=[]
for i in [15,18,17,19,16,20,21,22,23,24,25,26]:
 base=RAW/'EXERCISE C-2'/'ogma_long_exposures'
 f=base/f'CCD Image {i}.fit'
 if not f.exists():f=base/'iforgot'/f'CCD Image {i}.fit'
 if f.exists():
  im=fits.getdata(f);hd=fits.getheader(f);box=im[1338:1479,1842:1983];long.append(dict(image=i,t_s=hd['EXPTIME'],raw_max_DN=int(im.max()),stellar_region_peak_DN=int(box.max()),digital_clipped_pixels=int((im>=65535).sum())))
if long:export(pd.DataFrame(long),'ogma_raw_maxima.csv');summary['ogma_raw_maxima']=long
(OUT/'revision_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

# Export dark-region means without doing the spreadsheet calculations.
dark_rows=[]
for series,folder in [('time','Dark_Time'),('temperature','Dark_Temperature')]:
 for file in sorted((ASTRO/'A'/folder).glob('*.fit')):
  im,hd=fits.getdata(file,header=True)
  if series=='temperature' and (hd.get('EXPTIME')!=300 or hd.get('CCD-TEMP') is None):continue
  dark_rows.append({'series':series,'file':file.name,'raw_mean_DN':float(im[1478:1578,1478:1578].astype(float).mean()),'raw_SD_DN':float(im[1478:1578,1478:1578].astype(float).std(ddof=1)),'region_pixels':10000})
if dark_rows:export(pd.DataFrame(dark_rows),'dark_region_raw_means.csv')

# Repeatability image and profile used in the paper.
if (RAW/'EXERCISE C-2/Ogma-0001.fit').exists() and (RAW/'EXERCISE C-2/Ogma calibrated/Ogma-0001_calibratedogma.fits').exists():
 base=RAW/'EXERCISE C-2'
 raw=fits.getdata(base/'Ogma-0001.fit')
 # AIJ T1 is near (952, 828); its y coordinate is reflected relative to the FITS array.
 box=raw[660:740,910:990];by,bx=np.unravel_index(box.argmax(),box.shape);y,x=by+660,bx+910
 cut=raw[y-20:y+21,x-20:x+21]
 cal=fits.getdata(base/'Ogma calibrated/Ogma-0001_calibratedogma.fits')
 plt.rcParams.update({'font.size':11,'axes.titlesize':11,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10})
 fig,ax=plt.subplots(1,2,figsize=(6.5,3.0),layout='constrained')
 ax[0].imshow(cut,origin='upper',cmap='gray',extent=(-20.5,20.5,20.5,-20.5));ax[0].set(title='(a) Ogma core',xlabel='Offset (pixels)',ylabel='Offset (pixels)',xticks=[-20,0,20],yticks=[-20,0,20])
 star=cal[1207:1268,191:252];ax[1].imshow(star,origin='upper',cmap='gray',vmin=np.percentile(star,10),vmax=np.percentile(star,99.8));ax[1].plot(30,30,'+',color='#eead55',ms=7);ax[1].set(title='(b) Field star',xlabel='Cutout x (pixels)',ylabel='Cutout y (pixels)',xticks=[0,30,60],yticks=[0,30,60])
 fig.savefig(REPORT/'repeatability_evidence.pdf',bbox_inches='tight');plt.close(fig)
