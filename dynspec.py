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

import ismfit
import ism_estimates
import lmfit
import collections

class DynSpec:
    def __init__(self,source,epochs,telescope,tints,profile,onp,offp,freqs,on,off,guppi=False,alpha=-4.0,discard=0.01):
        self.guppi = guppi
        self.alpha = alpha
        self.source = source
        self.epochs = epochs
        self.tints = tints
        self.telescope = telescope
        self.profile = profile
        self.onp = onp
        self.freqs = freqs
        self.offp = offp
        mask = maskFreqs(self.freqs, self.telescope)
        self.freqs = self.freqs[mask]
        self.on = on[:,:,mask]
        self.off = off[:,:,mask]
        self.ds = cleanAndNormalize(self.on, self.off,discard=discard)
        if self.ds.shape[0] % 2 == 1:
            self.epochs = self.epochs[:-1]
            self.tints = self.tints[:-1]
            self.ds = self.ds[:-1,:]
        self.times = (self.epochs - self.epochs.min())*86400
        if self.times.shape[0] < 4:
            raise Exception("Not long enough")
        self.sds,self.fc,self.sfreqs = stretchDS(self.ds,self.freqs,alpha=alpha)
        self.acf = computeAcf(self.ds)
        self.sacf = computeAcf(self.sds)
        self.refit()
    def plot(self,fig=None,stretch=False):
        if stretch:
            plotDynSpecAcf(self.sds, self.sacf, self.sfreqs, self.times, df = (self.freqs[1]-self.freqs[0]),fig=fig, guppi=self.guppi, profile=self.profile, onp=self.onp, offp=self.offp,fit=self.fit_stretch)
        else:
            plotDynSpecAcf(self.ds, self.acf, self.freqs, self.times, fig=fig, guppi=self.guppi, profile=self.profile, onp=self.onp, offp=self.offp,fit = self.fit)
    def refit(self,tdif0=None,fdif0=None):
        if not(tdif0 or fdif0):
            print "looking up initial values"
            self.tdif0,self.fdif0 = ism_estimates.findEstimates(self.source, freq = self.freqs.mean(), bw = np.abs(self.freqs[0]-self.freqs[-1]))
            tdif0 = self.tdif0
            fdif0 = self.fdif0
        self.fit = AcfFit(dt = self.times[1]-self.times[0], df = self.freqs[1]-self.freqs[0],
                          acf= self.acf, tdif0 = tdif0, fdif0 = fdif0)
        self.fit_stretch = AcfFit(dt = self.times[1]-self.times[0], df = self.freqs[1]-self.freqs[0],acf = self.sacf)

def plotDynSpecAcf(ds,acf,frqs,times,df = None,fig = None,guppi=False,profile=None,onp=None,offp=None,fit=None):
    times = times - times[0]
    if df is None:
        df = np.abs(frqs[1]-frqs[0])
    dt = np.diff(times).mean()
    if guppi:
        if ds.shape[1] == 512:
            maxf = 100.0
        else:
            maxf = 50.0
    else:
        maxf = 3.0
    if fit is None:
        fit1 = AcfFit(dt,df,acf)
    else:
        fit1 = fit
    fit2 = AcfFit(dt,df,acf,maxf,fdif0=fit1.fdif0,tdif0=fit1.tdif0)
    if not guppi:
        fscale = int(5*fit2.params['fdif'].value/df)
        print fscale
        if fscale < 512:
            fscale = 512
    else:
        fscale = ds.shape[1]
    if fig is None:
        f = plt.figure()
    else:
        f = fig
    if profile is None:
        ax1 = f.add_subplot(3,2,1)
        ax2 = f.add_subplot(3,2,2)
    else:
        ax1 = f.add_subplot(3,3,1)
        ax2 = f.add_subplot(3,3,2)
        axp = f.add_subplot(3,3,3)
        nbin = profile.shape[0]
        bins = np.linspace(0,1,nbin)
        axp.plot(bins,profile,lw=2)
        if len(onp) == 2:
            axp.axvspan(onp[0]*1.0/nbin,onp[1]*1.0/nbin,color='g',alpha=0.4)
            axp.axvspan(offp[0]*1.0/nbin,offp[1]*1.0/nbin,color='r',alpha=0.4)
        else:
            axp.fill_between(bins, profile,where=onp,color='g')
            axp.fill_between(bins, profile,where=offp,color='r')
            axp.set_ylim(profile.min(),profile.max())
        for tl in axp.yaxis.get_ticklabels():
            tl.set_visible(False)
    ax3 = f.add_subplot(3,1,2)
    ax4 = f.add_subplot(3,1,3)
    ax1.plot(fit1.t,fit1.tacfn)
    ax1.plot(fit1.t,fit1.fittacfn, linewidth=2, alpha=0.5, label=((r'$\exp{(-(t/(%.2f\pm%.2f)^{5/3})}$' % (fit1.params['tdif'].value,fit1.params['tdif'].stderr))))
    ax1.legend(prop=dict(size='small'))
    ax1.set_title('Normalized temporal ACF')
    ax1.set_xlabel('Lag (s)')
    ax2.plot(fit1.fr,fit1.facfn)
    ax2.plot(fit1.fr,fit1.fitfacfn, linewidth=2, alpha=0.5, label=((r'$\exp{(\frac{-(\log{2})\Delta f}{%.4f\pm%.4f})}$' % (fit1.params['fdif'].value,fit1.params['fdif'].stderr))))
    if True:
        ax2.plot(fit2.fr,fit2.fitfacfn, linewidth=2, alpha=0.5, label=((r'$\exp{(\frac{-(\log{2})\Delta f}{%.4f\pm%.4f})}$' % (fit2.params['fdif'].value,fit2.params['fdif'].stderr))))
    if not guppi and (fscale*df < fit1.fr.max()):
        ax2.set_xlim(0,fscale*df)
    ax2.legend(prop=dict(size='small'))
    ax2.set_title('Normalized frequency ACF')
    ax2.set_xlabel('Lag (MHz)')
    if not guppi:
        ds2 = rebinAxis(ds, 2048, axis=1)
        ds3 = rebinAxis(ds2, 256, axis = 1)
        im = ax4.imshow(ds2,aspect='auto',origin='lower',extent=[frqs[0],frqs[-1],times[0],times[-1]])
        cs = ax4.contour(ds3,extent=[frqs[0],frqs[-1],times[0],times[-1]],
                 levels = ds3.max()*np.linspace(0.5,1,5), cmap=plt.cm.gray_r, linewidths=1,alpha=0.5)

    else:
        im = ax4.imshow(ds,aspect='auto',origin='lower',extent=[frqs[0],frqs[-1],times[0],times[-1]])
        cs = ax4.contour(ds,extent=[frqs[0],frqs[-1],times[0],times[-1]],
                 levels = ds.max()*np.linspace(0.2,1,5), cmap=plt.cm.gray_r, linewidths=1,alpha=0.5)
    ax4.text(0.1,0.9,"Dynamic Spectrum",transform=ax4.transAxes,bbox=dict(facecolor='w',alpha=0.5))
    ax4.set_xlabel('Freq (MHz)')
    ax4.set_ylabel('Time (s)')
    #im.set_clim(-0.01,0.05)
    cb = f.colorbar(im,ax=ax4)
    try:
        cb.add_lines(cs)
    except ValueError:
        pass
    
    if fscale > 1024:
        scaleby = int(np.ceil(fscale/1024))
        acf2 = rebinAxisBy(acf, scaleby, axis=1)
    else:
        acf2 = acf
        scaleby = 1
    print fscale,scaleby
    if not guppi:
        acf2 = fit1.normalize(acf2) #/acf2[tuple(np.array(acf2.shape)/2 + 1)]
    else:
        acf2 = fit1.normalize(acf2) #/acf2[tuple(np.array(acf2.shape)/2 + 1)]
#    acf2 = acf2 - acf2[int(0.75*acf2.shape[0]):,int(0.9*acf2.shape[1]):].mean()

#    acfs = acf2.flatten()
#    acfs.sort()
#    levels = [acfs[int(acfs.shape[0]*x)] for x in (1-10**-np.linspace(2,5,10))]
#    levels = np.linspace(0.4,1,6)
    nf = acf2.shape[1]
    print dt,df,acf2.shape
    if not guppi:
        extent = [-df*fscale, df*fscale, -dt*acf2.shape[0]/2.0,dt*acf2.shape[0]/2.0]
        im = ax3.imshow(acf2[:,nf/2-fscale/scaleby:nf/2+fscale/scaleby],aspect='auto',origin='lower',extent=extent)
    else:
        extent = [-df*nf/2.0, df*nf/2.0, -dt*acf2.shape[0]/2.0,dt*acf2.shape[0]/2.0]
        im = ax3.imshow(acf2,aspect='auto',origin='lower',extent=extent)
    if not guppi:
        smoothfactor =8
    else:
        smoothfactor =1
    acfsm = rebinAxisBy(acf2,smoothfactor,axis=1)
    levels = np.array([.5,.75,.9,.95,.99])
    levels = np.linspace(.1,1,10)
    nfsm = acfsm.shape[1]
    ax3.text(0.1,0.9,"2-D ACF",transform=ax3.transAxes,bbox=dict(facecolor='w',alpha=0.5))
    im.set_clim(-.1,1)
    fextent = fscale/(scaleby*smoothfactor)
    if not guppi:
        cs = ax3.contour(acfsm[:,nfsm/2-fextent:nfsm/2+fextent],extent=extent,
                         levels = levels, cmap=plt.cm.gray_r, linewidths=1,alpha=0.5)
    else:
        cs = ax3.contour(acfsm,extent=extent,
                         levels = levels, cmap=plt.cm.gray_r, linewidths=1,alpha=0.5)
    cb = f.colorbar(im,ax=ax3)
    try:
        cb.add_lines(cs)
    except ValueError:
        pass
    
    ax3.set_xlabel('Lag (MHz)')
    ax3.set_ylabel('Lag (s)')

def cleanAndNormalize(on,off,discard=0.01):
    #data dimensions are subint, poln, channel
    ds = on/off - 1.0
    ds[~np.isfinite(ds)]=0.0
    for pol in range(ds.shape[1]):
        pdata = ds[:,pol,:]
        shape = pdata.shape
        ndat = np.prod(shape)
        pdata = pdata.reshape((ndat,))
        ordering = pdata.argsort()
        if int(ndat*discard):
            pdata[ordering[:int(ndat*discard)]] = 0.0
            pdata[ordering[-int(ndat*discard):]] = 0.0
        pdata = pdata.reshape(shape)
        ds[:,pol,:] = pdata
        
    return ds.mean(1)  #scrunch the polns


def computeAcf(in1,correct=True):
    s1 = np.array(in1.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex))
    size = 2**np.ceil(np.log2(s1*2))
    IN1 = np.fft.fftn(in1,size)
    IN1 *= np.conj(IN1)
    ret = np.fft.ifftn(IN1)
    del IN1
    if not complex_result:
        ret = ret.real
    osize = s1
    output = scipy.signal.signaltools._centered(np.fft.fftshift(ret),osize)
    if not correct:
        return output
    corrections = []
    for d in range(output.ndim):
        ndat = output.shape[d]
        c = np.zeros((ndat,))
        half = ndat % 2
        c[:ndat/2] = (ndat*1.0/(ndat-np.arange(1,ndat/2+1)))[::-1]
        c[ndat/2:] = ndat*1.0/(ndat-np.arange(ndat/2 + half))
        corrections.append(c)
    return output * np.outer(corrections[0],corrections[1])

        
Param = collections.namedtuple('Param', ('name','value','stderr','correl','max','min'))
def makePickleableParam(param):
    return Param(*[getattr(param,x) for x in Param._fields])        
class AcfFit():
    def __init__(self,dt,df,acf,maxf=None,tdif0=None,fdif0=None):
        self.tdif0=tdif0
        self.fdif0=fdif0
        nt,nf = acf.shape
        self.facf = acf[nt/2+1,nf/2+1:]
        self.fr = df*np.arange(1,self.facf.shape[0]+1)
        if maxf:
            maxidx = int(maxf/df)-1
            self.fr = self.fr[:maxidx]
            self.facf = self.facf[:maxidx]
        self.tacf = acf[nt/2+1:,nf/2]
        self.t = dt*np.arange(1,self.tacf.shape[0]+1)
        if maxf is None:
            mi = ismfit.simFitAcf(self.t, self.tacf, self.fr, self.facf,tdif0=tdif0,fdif0=fdif0)
        else:
            mi = ismfit.fitIf3(self.fr, self.facf, fdif0=fdif0)
        try:
            self.ci,self.trace = lmfit.conf_interval(mi,trace=True)
        except:
            self.ci = None
            self.trace = None
        self.params = dict([(x,makePickleableParam(mi.params[x])) for x in mi.params.keys()])
        self.fitfacf = ismfit.gammaIf3(self.params,self.fr)
        self.offset = self.params['offs'].value
        self.scale = self.params['scale'].value
        if maxf is None:
            self.fittacf = ismfit.gammaIs3(self.params,self.t)
            self.fittacfn = self.normalize(self.fittacf)
        self.fitfacfn = self.normalize(self.fitfacf)
        self.facfn = self.normalize(self.facf)
        self.tacfn = self.normalize(self.tacf)
#        self.mi = mi
    def normalize(self,data):
        return (data-self.offset)/self.scale

def computeCcf(in1,in2):
    """Correlate two N-dimensional arrays using FFT. See convolve.

    """
    s1 = np.array(in1.shape)
    s2 = np.array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))
    size = s1+s2
    IN1 = np.fft.fftn(in1,size)
    IN1 *= np.conj(np.fft.fftn(in2,size))
    ret = np.fft.ifftn(IN1)
    del IN1
    if not complex_result:
        ret = ret.real
    if np.product(s1,axis=0) > np.product(s2,axis=0):
        osize = s1
    else:
        osize = s2
    return scipy.signal.signaltools._centered(ret,osize)

def rebinAxisBy(data,nfact,axis=0,oper=np.mean):
    nout = int(data.shape[axis]/nfact)
    return rebinAxis(data,nout,axis=axis,oper=oper)
    
def rebinAxis(data,nout,axis=0, oper = np.mean):
    nout = int(nout)
    nstart = data.shape[axis]
    nchk = nstart/nout
    nstop = nchk*nout
    index = [slice(None) for x in range(axis)] + [slice(nstop)] 
    if data.ndim > axis:
        index = index + [slice(None) for x in range(data.ndim-axis-1)]
    index = tuple(index)
    newshape = list(data.shape)
    newshape[axis] = nout
    newshape.insert(axis+1,-1)
    return np.mean(data[index].reshape(tuple(newshape)),axis=axis+1)

def stretchDS(ds,f,alpha=-4):
    fc = np.sqrt(f.max()*f.min())
    f2 = np.cumsum((f/fc)**alpha)
    nfout = np.floor(f2.max())
    rds = np.zeros((ds.shape[0],nfout))
    x = np.arange(nfout)
    for k in range(ds.shape[0]):
        rds[k,:] = np.interp(x,f2,ds[k,:])
    fout = np.interp(x,f2,f)
    return rds,fc,fout

def findOnOff(prof):
    med = np.median(prof)
    offrms = prof[prof<med].std()
    onreg = (prof>(med + 5*offrms)).astype('int')
    onreg = (prof>(prof.min()+0.15*prof.ptp()))
    offreg = ~onreg
    return onreg,offreg

def maskFreqs(freqs,telescope):
    if telescope == 'Arecibo':
        mask = (((freqs > 1100) & (freqs < 1790)) |
                ((freqs > 1720) & (freqs < 2410)) |
                ((freqs > 420) & (freqs < 445)) |
                ((freqs > 302) & (freqs < 350)))
        return mask
    else:
        mask = (((freqs > 1150) & (freqs< 1880)) | (freqs <920))
        return mask 
    
def pickle(fn,obj):
    fh = open(fn,'w')
    cPickle.dump(obj, fh, protocol=cPickle.HIGHEST_PROTOCOL)
    fh.close()
    
def unpickle(fn):
    fh = open(fn,'r')
    res = cPickle.load(fh)
    fh.close()
    return res