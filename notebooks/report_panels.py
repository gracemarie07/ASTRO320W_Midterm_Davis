"""Compact report panels from recorded measurements; does not remeasure the observations."""
from pathlib import Path
import os, shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data';R=ROOT/'report';O=ROOT/'results'
RAW=Path(os.environ.get('ASTRO320W_DATA',str(ROOT/'raw')))
plt.rcParams.update({'font.size':11.5,'axes.labelsize':11.5,'axes.titlesize':11.5,'legend.fontsize':11,'xtick.labelsize':11,'ytick.labelsize':11})
def save(f,n):
 f.savefig(O/(n+'.pdf'),bbox_inches='tight');f.savefig(O/(n+'.png'),dpi=200,bbox_inches='tight');plt.close(f);shutil.copy2(O/(n+'.pdf'),R/(n+'.pdf'))
# Three complementary views of the dark data.
a=pd.read_csv(D/'dark_time_with_errors.csv');b=pd.read_csv(D/'dark_temperature_with_errors.csv')
f=plt.figure(figsize=(6.5,6.1),layout='constrained');gs=f.add_gridspec(2,2);axs=[f.add_subplot(gs[0,:]),f.add_subplot(gs[1,0]),f.add_subplot(gs[1,1])]
x=a.t_s.to_numpy();q=a.Q_e.to_numpy();coef=np.polyfit(x[1:],q[1:],1);xx=np.linspace(0,62,100)
axs[0].errorbar(x[1:],q[1:],yerr=a.Q_formal_error_e.iloc[1:],fmt='o');axs[0].plot(xx,np.polyval(coef,xx));axs[0].plot(x[0],q[0],'s',mfc='none',label='0.1 s: excluded');axs[0].set(xlabel='Exposure time (s)',ylabel='Charge (electrons/pixel)',title='(a) Exposure-time dependence');axs[0].legend()
axs[1].errorbar(b.T_start_C,b.Q_e,yerr=b.Q_formal_error_e,fmt='o');axs[1].set(xlabel='Start temperature (C)',ylabel='Charge (electrons/pixel)',title='(b) Temperature series')
x=1000/(b.T_start_C+273.15);y=np.log(b.I_e_s);coef=np.polyfit(x,y,1);axs[2].plot(x,y,'o');xx=np.linspace(x.min(),x.max(),100);axs[2].plot(xx,np.polyval(coef,xx));axs[2].set(xlabel='1000 / T (1/K)',ylabel='ln(I / Iref)',title='(c) Thermal fit')
for ax in axs:ax.grid(alpha=.2)
save(f,'dark_panels')
# Paired PTCs; use the saved numeric sources explicitly.
a=pd.read_csv(D/'original_aij_sky_ptc.csv');b=pd.read_csv(D/'dome_ptc_digitized.csv')
f,axs=plt.subplots(1,2,figsize=(6.5,3.1),layout='constrained')
for ax,x,y,title in [(axs[0],a.mean_DN,a.spatial_SD_DN,'(a) Sky flats'),(axs[1],b.mean_DN_approx,b.spatial_SD_DN_approx,'(b) Dome flats (digitized)')]:
 ax.loglog(x,y,'o',ms=4);ax.set(xlabel='Mean signal (DN)',ylabel='Spatial SD (DN)',title=title);ax.grid(alpha=.2,which='both');ax.tick_params(axis='both',which='minor',labelsize=11)
save(f,'ptc_panels')
# Focus plus observed raw profile.
f,axs=plt.subplots(1,2,figsize=(6.5,3.2),layout='constrained')
x=np.array([17000,17000,17000,16500,16000,16000,15500,15000]);y=np.array([11.55,12.41,11.19,10.30,9.75,10.68,10.72,12.74]);axs[0].plot(x,y,'o',ms=4)
c=np.polyfit(x-16000,y,2);xx=np.linspace(15000,17000,100);axs[0].plot(xx,np.polyval(c,xx-16000),'--',lw=1)
for p in np.unique(x):
 v=y[x==p]
 if len(v)>1:axs[0].errorbar(p,v.mean(),yerr=v.std(ddof=1),fmt='s',color='#ad622b',capsize=3)
axs[0].set(xlabel='Focus position (micrometers)',ylabel='FWHM (binned pixels)',title='(a) Focus sequence');axs[0].tick_params(axis='x',labelrotation=25)
p=RAW/'EXERCISE C-2'/'CCD Image 6.fit'
if not p.exists():
 raise FileNotFoundError('Set ASTRO320W_DATA to the raw collection containing EXERCISE C-2/CCD Image 6.fit')
if p.exists():
 im=fits.getdata(p);iy,ix=np.unravel_index(im.argmax(),im.shape);axs[1].plot(np.arange(-20,21),im[iy,ix-20:ix+21]);axs[1].set(xlabel='Offset from peak (pixels)',ylabel='Raw signal (DN)',title='(b) Raw focus Image 6')
for ax in axs:ax.grid(alpha=.2)
save(f,'focus_panels')
# Repeatability panels show the observed scatter, not nominal errors alone.
a=pd.read_csv(D/'ogma.csv').iloc[:15];rates=a['Source-Sky_T1']/a.EXPTIME;b=pd.read_csv(D/'unsaturated_field_star.csv');r=b.Rate_DN_s
f,axs=plt.subplots(1,2,figsize=(6.5,3.1),layout='constrained')
for ax,v,title in [(axs[0],rates,'(a) Ogma: nonlinear cores'),(axs[1],r,'(b) Unsaturated field star')]:
 dev=100*(v-v.mean())/v.mean();sc=dev.std(ddof=1);ax.plot(np.arange(1,16),dev,'o-',ms=3);ax.axhline(0,color='.4',lw=.7);ax.axhspan(-sc,sc,color='#245884',alpha=.1);ax.set(xlabel='Exposure number',ylabel='Deviation from mean (%)',title=title,xticks=[1,5,10,15]);ax.grid(alpha=.2)
save(f,'repeatability_panels')
# Exposure-time table values are the original AIJ 30-pixel aperture measurements.
t=np.array([6,8,10,12,15,20,45,60,90,120,150,300]);w=[6.79,6.64,7.30,7.83,8.61,9.22,12.50,14.50,16.36,17.50,17.24,22.97];rate=np.array([297312,299378,275574,242015,237392,197273,159497,158589,131787,113424,87429,69921])
f,axs=plt.subplots(1,2,figsize=(6.5,3.1),layout='constrained');axs[0].plot(t,w,'o-',ms=3);axs[0].set(ylabel='FWHM (unbinned pixels)',title='(a) Image width');axs[1].plot(t,rate/1e5,'o-',ms=3);axs[1].set(ylabel='Count rate (100000 DN/s)',title='(b) Count rate')
for ax in axs:ax.set_xlabel('Exposure time (s)');ax.grid(alpha=.2)
save(f,'exposure_panels')
# Two color diagrams: full range and central enlargement with formal covariance.
c=pd.read_csv(D/'m29_reconciled.csv');c=c[c.included]
f,axs=plt.subplots(1,2,figsize=(6.5,3.4),layout='constrained')
for ax,title in zip(axs,['(a) All ten stars','(b) Central group']):
 for _,r in c.iterrows():
  if ax is axs[1] and r.Star==8:continue
  ax.errorbar(r.BV,r.VR,xerr=r.BV_error,yerr=r.VR_error,fmt='o',capsize=2,ms=3,color='#245884')
  offsets={1:(-6,-15),7:(-11,5),3:(6,-13),4:(-13,5),5:(-8,-16),12:(6,-12),2:(6,3),6:(-12,-14),9:(-11,5),8:(-12,-14)}
  ax.annotate(str(int(r.Star)),(r.BV,r.VR),xytext=offsets[r.Star],textcoords='offset points',fontsize=11)
  if ax is axs[1]:
   cv=np.array([[r.BV_error**2,r.covariance],[r.covariance,r.VR_error**2]]);v,u=np.linalg.eigh(cv);ang=np.degrees(np.arctan2(u[1,1],u[0,1]));ax.add_patch(Ellipse((r.BV,r.VR),2*np.sqrt(v[1]),2*np.sqrt(v[0]),angle=ang,facecolor='#245884',alpha=.25))
 ax.set(xlabel='Instrumental B - V (mag)',ylabel='Instrumental V - R (mag)',title=title);ax.grid(alpha=.2)
axs[0].set(xlim=(.76,1.60),ylim=(.15,.63));axs[1].set(xlim=(.79,1.12),ylim=(.17,.405));save(f,'m29_color_panels')
print('Generated six compact report panel figures.')

# Identification figures at their final report dimensions.
if os.environ.get('ASTRO320W_ORGANIZED') and (RAW/'EXERCISE C-2/Ogma calibrated').exists():
 root=ROOT;r=root/'report';d=pd.read_csv(root/'data/m29_reconciled.csv')
 im=fits.getdata(Path(os.environ['ASTRO320W_ORGANIZED'])/'D/AIJ_Reference_Stacks/M29_V_stack_bias_corrected.fit');lo,hi=np.percentile(im,[15,99.8]);f,ax=plt.subplots(figsize=(4.8,4.6),layout='constrained')
 ax.imshow(np.arcsinh(np.maximum(im-lo,0)/((hi-lo)/12)),cmap='gray',origin='upper',vmin=0,vmax=np.arcsinh(12))
 for _,v in d.iterrows():
  y=1528-v.y_IJ;clr='#73efff' if v.included else '#ffab40';ax.add_patch(plt.Circle((v.x,y),16,fill=False,color=clr,lw=1));ax.annotate(str(int(v.Star)),(v.x,y),xytext=(9,-3),textcoords='offset points',color=clr,fontsize=12,weight='bold')
 ax.set(xlabel='Array x (pixels)',ylabel='Array y (pixels)');ax.tick_params(labelsize=11);ax.xaxis.label.set_size(12);ax.yaxis.label.set_size(12);f.savefig(r/'M29_finding_chart.pdf',bbox_inches='tight');plt.close(f)
 files=sorted((RAW/'EXERCISE C-2/Ogma calibrated').glob('*.fits'));im=fits.getdata(files[0]);cut=im[1207:1268,191:252];f,ax=plt.subplots(figsize=(3.1,2.6),layout='constrained');ax.imshow(cut,cmap='gray',origin='upper',vmin=np.percentile(cut,10),vmax=np.percentile(cut,99.8));ax.plot(30,30,'+',color='#eead55',ms=7);ax.set(xlabel='Cutout x (pixels)',ylabel='Cutout y (pixels)');ax.tick_params(labelsize=11);ax.xaxis.label.set_size(12);ax.yaxis.label.set_size(12);f.savefig(r/'unsaturated_field_star_identification.pdf',bbox_inches='tight');plt.close(f)
