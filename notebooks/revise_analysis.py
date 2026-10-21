"""Reproduce the report revisions. Raw observations are optional; see README."""
from pathlib import Path
import os, json, csv, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import linregress
import openpyxl
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
# M29: filter labels determine the order; identifiers retain the matching AIJ export.
df=pd.read_csv(DATA/'M29.csv')
frames={b:df[df.Label.str.contains('M29_'+b+'_')].iloc[0] for b in 'BVR'}
rows=[]
for star in range(1,13):
 suffix='T1' if star==1 else f'C{star}'
 vals={b:float(frames[b]['Source-Sky_'+suffix]) for b in 'BVR'}
 errs={b:float(frames[b]['Source_Error_'+suffix]) for b in 'BVR'}
 k=2.5/np.log(10);bv=-2.5*np.log10(vals['B']/vals['V']);vr=-2.5*np.log10(vals['V']/vals['R'])
 rows.append(dict(Star=star,AIJ_ID=suffix,**{b:vals[b] for b in 'BVR'},**{b+'_error':errs[b] for b in 'BVR'},BV=bv,VR=vr,BV_error=k*np.hypot(errs['B']/vals['B'],errs['V']/vals['V']),VR_error=k*np.hypot(errs['V']/vals['V'],errs['R']/vals['R']),covariance=-(k*errs['V']/vals['V'])**2,R_peak=float(frames['R']['Peak_'+suffix]),x=float(frames['V']['X(IJ)_'+suffix]),y_IJ=float(frames['V']['Y(IJ)_'+suffix]),included=star not in [10,11]))
colors=pd.DataFrame(rows);export(colors,'m29_reconciled.csv')
c=colors[colors.included]
fig,ax=plt.subplots(figsize=(6.5,4.2),layout='constrained')
ax.errorbar(c.BV,c.VR,xerr=c.BV_error,yerr=c.VR_error,fmt='o',capsize=3,color='#245884')
for _,r in c.iterrows():ax.annotate(str(int(r.Star)),(r.BV,r.VR),xytext=(5,7 if r.Star not in [3,5,12] else -13),textcoords='offset points',fontsize=10)
ax.set(xlabel='Instrumental B - V (mag)',ylabel='Instrumental V - R (mag)',xlim=(.76,1.60),ylim=(.15,.63));ax.grid(alpha=.2);save(fig,'m29_color_color')
# Second view explicitly reveals the shared-V measurement covariance.
from matplotlib.patches import Ellipse
fig,ax=plt.subplots(figsize=(6.5,4.0),layout='constrained')
for _,r in c[c.Star!=8].iterrows():
 cov=np.array([[r.BV_error**2,r.covariance],[r.covariance,r.VR_error**2]])
 v,u=np.linalg.eigh(cov);ang=np.degrees(np.arctan2(u[1,1],u[0,1]))
 ax.add_patch(Ellipse((r.BV,r.VR),2*np.sqrt(v[1]),2*np.sqrt(v[0]),angle=ang,facecolor='#245884',alpha=.25,edgecolor='#245884'))
 ax.plot(r.BV,r.VR,'.',color='#245884');ax.annotate(str(int(r.Star)),(r.BV,r.VR),xytext=(5,7 if r.Star not in [3,5,12] else -13),textcoords='offset points',fontsize=10)
ax.set(xlabel='Instrumental B - V (mag)',ylabel='Instrumental V - R (mag)',xlim=(.79,1.12),ylim=(.17,.40));ax.grid(alpha=.2);save(fig,'m29_color_detail')
stack=ASTRO/'D/AIJ_Reference_Stacks/M29_V_stack_bias_corrected.fit'
if stack.exists():
 im=fits.getdata(stack).astype(float);fig,ax=plt.subplots(figsize=(6.5,6),layout='constrained')
 lo,hi=np.percentile(im,[15,99.8]);ax.imshow(np.arcsinh(np.maximum(im-lo,0)/max((hi-lo)/12,1)),origin='upper',cmap='gray',vmin=0,vmax=np.arcsinh(12))
 for _,r in colors.iterrows():
  y=im.shape[0]-r.y_IJ;clr='#ffab40' if not r.included else '#73efff'
  ax.add_patch(plt.Circle((r.x,y),15,fill=False,edgecolor=clr,lw=1.2));ax.annotate(str(int(r.Star)),(r.x,y),xytext=(17,-3),textcoords='offset points',color=clr,fontsize=11,weight='bold')
 ax.set(xlabel='Array x (pixels)',ylabel='Array y (pixels)');save(fig,'M29_finding_chart')
# Joint paired-flat gain regression, including the low-signal read-noise term.
faint=pd.read_csv(DATA/'dome_0p2s_gain_pairs_1792460302773.csv')
wb=openpyxl.load_workbook(next((ROOT/'data').glob('*.xlsx')),data_only=True)
bright=[r[:4] for r in wb.active.iter_rows(min_row=2,values_only=True) if r[0] is not None and all(isinstance(x,(int,float)) for x in r[1:4])]
S=list((faint.F1_mean_DN+faint.F2_mean_DN)/2)+[(r[1]+r[2])/2 for r in bright]
V=list(faint.difference_variance_DN2/2)+[r[3]**2/2 for r in bright]
fit=linregress(S,V);summary['gain']={'value':1/fit.slope,'sem':fit.stderr/fit.slope**2,'intercept_DN2':fit.intercept}
export(pd.DataFrame({'S_DN':S,'V_DN2':V}),'gain_regression_inputs.csv')
fig,ax=plt.subplots(figsize=(6.5,4.1),layout='constrained');ax.plot(S,V,'o',label='Disjoint flat pairs');xx=np.linspace(0,max(S)*1.03,100);ax.plot(xx,fit.intercept+fit.slope*xx,label='Variance = signal/gain + intercept');ax.set(xlabel='Mean paired signal (DN)',ylabel='Half difference-image variance (DN squared)');ax.legend();ax.grid(alpha=.2);save(fig,'gain_regression')
# Plate scale and full observing-status record.
e=pd.read_csv(DATA/'E_verified_repeats.csv');means=[];fig,ax=plt.subplots(figsize=(6.5,4.3),layout='constrained')
offsets={'HR6200':(-15,-23),'HR8815':(10,7),'HR7310':(9,-5),'HR7751':(-65,8)}
for target,g in e.groupby('target',sort=False):
 d=np.hypot(g.x2-g.x1,g.y2-g.y1);theta=g.catalog_arcsec.iloc[0];scales=theta/d;means.append(scales.mean())
 ax.errorbar(theta,d.mean(),yerr=d.std(ddof=1)/np.sqrt(len(d)),fmt='o',capsize=3);ax.annotate(target,(theta,d.mean()),xytext=offsets[target],textcoords='offset points',fontsize=10)
p=np.mean(means);se=np.std(means,ddof=1)/2;xx=np.linspace(0,225,100);ax.plot(xx,xx/p,'--',color='.4',label='Equal-binary mean scale');ax.set(xlabel='Catalog separation (arcseconds)',ylabel='Measured separation (binned pixels)',xlim=(0,225),ylim=(0,185));ax.legend();ax.grid(alpha=.2);save(fig,'E_separation_pixels_vs_arcsec');summary['plate_scale']={'value':p,'sem':se,'FOV_arcmin':p*1528/60}
# Full-detector noise maps: disjoint-pair noise and mean single-frame spatial SD.
files=sorted(BIAS.glob('biasframe-*.fit'))
if len(files)==20:
 stack=np.stack([fits.getdata(f).astype(np.float32) for f in files]);master=stack.mean(axis=0);h,w=master.shape;grid=[];RN=np.zeros((8,8));SP=np.zeros((8,8))
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
 # Formal uncertainty in the central mean, with matched-noise assumptions stated in report.
 B=float(master[1478:1578,1478:1578].mean());rn=10.045;n=10000
 for subfolder,prefix in [('Dark_Time','dark_time_with_errors'),('Dark_Temperature','dark_temperature_with_errors')]:
  out=[]
  for f in sorted((ASTRO/'A'/subfolder).glob('*.fit')):
   hd=fits.getheader(f);t=hd.get('EXPTIME',0)
   if subfolder=='Dark_Temperature' and (t!=300 or hd.get('CCD-TEMP') is None):continue
   arr=fits.getdata(f)[1478:1578,1478:1578].astype(float);net=float(arr.mean()-B);q=1.34*net
   # Pixel-independence model: source Poisson + read noise + finite master bias.
   sd=1.34*np.sqrt((max(net,0)/1.34+rn**2*(1+1/20))/n)
   out.append(dict(file=f.name,t_s=t,T_start_C=hd.get('CCD-TEMP'),net_DN=net,Q_e=q,I_e_s=q/t,Q_formal_error_e=sd,I_formal_error_e_s=sd/t))
  export(pd.DataFrame(out).sort_values('t_s' if subfolder=='Dark_Time' else 'T_start_C'),prefix+'.csv')
 summary['central_bias_mean_DN']=B
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
