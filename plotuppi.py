import numpy as np
import glob
import os
import scipy.signal
import scipy.optimize
import matplotlib
from matplotlib import pyplot as plt
#matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 
import psrchive
import cPickle
import shutil
import dynspec

matplotlib.rcParams['xtick.labelsize']='small'
matplotlib.rcParams['ytick.labelsize']='small'

gbtpsrs = ['B1937+21',
           'J0218+4232',
           'J0340+4130',
           'J0613-0200',
           'J0645+5158',
#           'J0931-1902',
           'J1012+5307',
           'J1024-0719',
#           'J1327-0755',
#           'J1400-1438',
           'J1455-3330',
           'J1600-3053',
           'J1614-2230',
           'J1643-1224',
           'J1713+0747',
           'J1744-1134',
           'J1747-4036',
           'J1832-0836',
           'J1909-3744',
           'J1918-0642',
           'J2010-1323',
           'J2145-0750',
           'J2302+4442']

aopsrs = ['0023+0923',
 '0030+0451',
 '1640+2224',
 '1713+0747',
 '1738+0333',
 '1741+1351',
 '1853+1303',
 '1855+09',
 '1903+0327',
 '1910+1256',
 '1923+2515',
 '1937+21',
 '1944+0907',
 '1949+3106',
 '1953+29',
 '1955+2527',
 '2017+0603',
 '2043+1711',
 '2214+3000',
 '2317+1439']

def plotPulsar(psr, maxn = 6, band=820):
    plotdir = '/home/gjones/plots/aopanels'
    files = glob.glob('/lakitu/scratch/nanograv/uppi/analysis/dynspec/pkls/ds_*%s_*_ao_puppi.pkl' % psr)
    files.sort()
    if not files:
        print "no files found!"
        return
    dsfig = Figure(figsize=(8,8))
    dsfig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
    acffig = Figure(figsize=(8,8))
    acffig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
    ssfig = Figure(figsize=(8,8))
    ssfig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
    plotPulsarFiles(files,maxn=maxn,band=band,dsfig=dsfig,acffig=acffig,ssfig=ssfig)
    canvas = FigureCanvasAgg(dsfig)
    canvas.print_figure(os.path.join(plotdir,('%s_ds_%d.png' % (psr,band))))
    canvas = FigureCanvasAgg(acffig)
    canvas.print_figure(os.path.join(plotdir,('%s_acf_%d.png'% (psr,band))))
    canvas = FigureCanvasAgg(ssfig)
    canvas.print_figure(os.path.join(plotdir,('%s_ss_%d.png'% (psr,band))))
    
    
def plotSummary(psr,band=820,max = 300, maxerr=50):
    plotdir = '/home/gjones/plots/aopanels'
    files = glob.glob('/lakitu/scratch/nanograv/uppi/analysis/dynspec/pkls/ds_*%s_*_ao_puppi.pkl' % psr)
    files.sort()
    if not files:
        print "no files found!"
        return
    dss = [dynspec.unpickle(fn) for fn in files]
    dss = [ds for ds in dss if np.abs(ds.fc-band) < 100]
    fdif = np.array([ds.fit.params['fdif'].value for ds in dss])
    fdiferr = np.array([ds.fit.params['fdif'].stderr for ds in dss])
    epochs = np.array([ds.epochs.mean() for ds in dss])
    bad = (fdif > max) | (fdiferr > maxerr)
    epochs = epochs[~bad]
    fdif = fdif[~bad]
    fdiferr = fdiferr[~bad]
    fig = Figure(figsize=(6,3))
    fig.subplots_adjust(left=0.09,bottom=0.13,top=0.9,right=0.95)
    ax = fig.add_subplot(111)
    try:
        ax.errorbar(epochs,fdif,yerr=fdiferr,lw=0.5,mew=2,fmt='o',linestyle='-',elinewidth=1.5)
    except Exception, e:
        print e
    ax.grid(True)
    ax.set_xlabel('MJD')
    ax.set_ylabel('MHz')
    ax.set_ylim(0,5*np.median(fdif))#1.05*max)
    ds = dss[0]
    fig.suptitle("%s @ %s %d MHz Diffractive Bandwidth Estimates" % (psr,ds.telescope,band))
    canvas = FigureCanvasAgg(fig)
    canvas.print_figure(os.path.join(plotdir,('%s_fdif_%d.pdf'% (psr,band))))
    canvas.print_figure(os.path.join(plotdir,('%s_fdif_%d.png'% (psr,band))))
    
def plotPulsarFiles(files,maxn=6,band=820,dsfig=None,acffig=None,ssfig=None,discard=1e-2):
    if dsfig is None:
        dsfig = plt.figure()
    if acffig is None:
        acffig = plt.figure()
        
    axnum = 1
    maxax = maxn**2
    psr = None
    for fn in files:
        print fn
        ds = dynspec.unpickle(fn)
        if np.abs(ds.fc-band) > 100:
            print "skipping, fc is %.2f MHz" % ds.fc
            continue
        if psr is None:
            psr = ds.source
            telescope=ds.telescope
        ds2 = dynspec.cleanAndNormalize(ds.on,ds.off,discard=discard)
        ax = dsfig.add_subplot(maxn,maxn,axnum)
        ax.imshow(ds2,aspect='auto',origin='lower',extent=[ds.freqs[0],ds.freqs[-1],ds.times[0]/60.0,ds.times[-1]/60.0])
        ax.text(0.99,0.95,('%.0f' % ds.epochs[0]),ha='right',va='top',fontdict=dict(size='x-small'),transform=ax.transAxes,bbox=dict(alpha=0.5,fc='w'))
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax = acffig.add_subplot(maxn,maxn,axnum)
        ax.imshow(ds.acf,aspect='auto',origin='lower')
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.text(0.99,0.95,('%.0f' % ds.epochs[0]),ha='right',va='top',fontdict=dict(size='x-small'),transform=ax.transAxes,bbox=dict(alpha=0.5,fc='w'))
        
        ss = np.abs(np.fft.fftshift(np.fft.fft2(ds2)))
        ss = ss-ss[0].mean()
        ss = ss/ss.max()
        
        ax = ssfig.add_subplot(maxn,maxn,axnum)
        im = ax.imshow(20*np.log10(ss),aspect='auto',origin='lower',cmap=plt.cm.gray_r)
        im.set_clim(-60,0)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.text(0.99,0.95,('%.0f' % ds.epochs[0]),ha='right',va='top',fontdict=dict(size='x-small'),transform=ax.transAxes,bbox=dict(alpha=0.5,fc='w'))
        
        axnum+=1
        if axnum > maxax:
            break
    dsfig.suptitle("%s @ %s %d MHz Dynamic Spectra" % (psr,telescope,band))
    acffig.suptitle("%s @ %s %d MHz ACFs" % (psr,telescope,band))
    ssfig.suptitle("%s @ %s %d MHz Secondary Spectra" % (psr,telescope,band))
    