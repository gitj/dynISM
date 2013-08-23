"""
Routines for loading *uppi data files
"""

import os
import glob
import re

import numpy as np

from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 

import psrchive
import re

import dynspec

def processOne(fn,done,telescope,outdir,plotdir,redo=False):
    try:
        info = decodeFilename(os.path.split(fn)[1])
    except:
        print "couldn't decode",fn
        return fn,"couldn't decode"
    print info
    pklname = ('ds_%s_%s_%s_%s_%s.pkl' % (info['source'],info['mjd'],
                                      info['scan'],telescope,info['instrument']))
    outfn = os.path.join(outdir,'pkls',pklname)
    if (not redo) and (outfn in done):
        return fn,None
    fns = glob.glob(fn[:-9]+'*')
    
    print fns
    try:
#    if True:
        ds = loadDynSpecFromGPuppi(fns)
        print "loaded dynspec"
        dynspec.pickle(outfn, ds)
        print "pickled",outfn
        fig = Figure(figsize=(10,12))
        fig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
        ds.plot(fig=fig)
        plotname = os.path.join(plotdir,pklname + '.png')
        esc_fname = outfn.replace('_',r'\_')
        fig.suptitle(('%s @ %s %s' % (info['source'],telescope,esc_fname)),size='medium')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(plotname)

        fig = Figure(figsize=(10,12))
        fig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
        ds.plot(fig=fig,stretch=True)
        plotname = os.path.join(plotdir,'stretch',pklname + '.png')
        esc_fname = outfn.replace('_',r'\_')
        fig.suptitle(('%s @ %s %s' % (info['source'],telescope,esc_fname)),size='medium')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(plotname)
    except Exception,e:
        print fn,e
        return fn,e
    return fn,None

fndecode = re.compile(r"(?P<instrument>.uppi)_(?P<mjd>\d*)_(?P<source>[BJ]*\d{4}[+-]\d{2,4})_(?P<scan>\d*)_(?P<fileno>\d*)\.fits")

def decodeFilename(fn):
    res = fndecode.search(fn)
    res = res.groupdict()
    return res


def loadDynSpecFromGPuppi(fns,thresh=0.2,alpha=-4.0,discard=0.001):
    fns.sort()
    ars = [psrchive.Archive_load(fn) for fn in fns]
    source = ars[0].get_source()
    epochs = []
    times = []
    tints = []
    for ar in ars:
        nsubint = ar.get_nsubint()
        for x in range(nsubint):
            subint = ar.get_Integration(x)
            epochs.append(subint.get_start_time().in_days())
            times.append(subint.get_epoch().in_seconds())
            tints.append(subint.get_duration())
    times = np.array(times)
    tints = np.array(tints)
    epochs = np.array(epochs)
    times -= times[0]
    if times.shape[0] < 4:
        raise Exception("Observation is not long enough %s" % fn)
    telescope = ars[0].get_telescope()
    d = []
    weights = []
    for ar in ars:
        ar.bscrunch_to_nbin(256)
        ar.dedisperse()
        d.append(ar.get_data()*(ar.get_weights()[:,None,:,None]+1e-6))
    d = np.concatenate(d,axis=0)
    profile = d.mean(0)[0].mean(0)
    onp = profile > np.median(profile) + (np.max(profile)-np.median(profile))*thresh
    offp = profile < np.median(profile) + (np.max(profile)-np.median(profile))*0.1
    i0 = ars[0].get_Integration(0)
    freqs = np.array([i0.get_Profile(0,x).get_centre_frequency() for x in range(d.shape[2])])
    mask = dynspec.maskFreqs(freqs, telescope)
    freqs = freqs[mask]
    d = d[:,:,mask,:]
    if freqs[0] > freqs[1]:
        freqs = freqs[::-1]
        d = d[:,:,::-1,:]
    #self.on = d[:,:2,:,self.onp[0]:self.onp[1]].mean(3)
    #self.off = d[:,:2,:,self.offp[0]:self.offp[1]].mean(3)
    on = d[:,:2,:,onp].mean(3)
    off = d[:,:2,:,offp].mean(3)
    ds = dynspec.cleanAndNormalize(on, off,discard=discard)
    if ds.shape[0] % 2 == 1:
        times = times[:-1]
        ds = ds[:-1,:]
            
    ds = dynspec.DynSpec(source=source,
                         epochs=epochs,
                         telescope=telescope,
                         tints=tints,
                         profile=profile,
                         onp=onp,
                         offp=offp,
                         freqs=freqs,
                         on=on,
                         off=off,
                         guppi=True,alpha=alpha,discard=discard)
    return ds
            