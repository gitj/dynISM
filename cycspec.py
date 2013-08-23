"""
Routines for loading cycSpec data files
"""

import os
import glob
import re

import numpy as np

from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 

import psrchive

import dynspec

GPU_BASEDIRS = {'gbt':'/lakitu/data/scratch/nanograv/cs/gbt',
                'ao':'/lakitu/data/P2721'}

FILE_PATTERNS = {'gbt':'/lakitu/data/scratch/nanograv/cs/gbt/gpu%02d/rtcs_%d_*_%04d_*nth4_*.ar.fix',
                 'ao':'/lakitu/data/P2721/gpu%02d/rt/rtcs_%d_*_%04d_*nth4_*.ar.fix',
                 }

fndecode = re.compile(r"rtcs_(?P<mjd>\d*)_(?P<source>[BJ]*\d{4}[+-]\d{2,4})_(?P<scan>\d*)_.*_nth4_(?P<fileno>\d*)\.ar\.fix")

def processOne(fn,done,gpus,telescope,outdir,plotdir,redo=False):
    print fn
    info = decodeFilename(fn)
    print info
    pklname = ('ds_%s_%s_%s_%s.pkl' % (info['source'],info['mjd'],
                                      info['scan'],telescope))
    outfn = os.path.join(outdir,'pkls',pklname)
    if (not redo) and (outfn in done):
        return fn,None
    try:
    #if True:
        ds = loadDynSpecFromCycSpecScan(int(info['scan']), int(info['mjd']), 
                                        gpus=gpus, telescope=telescope)
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
    except Exception,e:
        print fn,e
        return fn,e
    return fn,None
                             

def decodeFilename(fn):
    res = fndecode.search(fn)
    res = res.groupdict()
    return res
def getProfileForScan(scanno,mjd,telescope='gbt'):
    files = getFilesForScan(scanno, mjd,telescope)
    return getProfile(files)
def getProfile(files):
    pfs = []
    for gpu in files.keys():
        for fn in files[gpu]:
            print "loading profile from:",fn
            ar = psrchive.Archive_load(fn)
#            ar.bscrunch_to_nbin(64)
            ar.tscrunch_to_nsub(1)
            ar.fscrunch_to_nchan(1)
            d = ar.get_data().squeeze().mean(0)
            pfs.append(d)
    return np.array(pfs).mean(0)

def getFreqMap(files,npfb=32):
    gpus = files.keys()
    freqs = {}
    for gpu in gpus:
        ar = psrchive.Archive_load(files[gpu][0])
        i0 = ar.get_Integration(0)
        nch = ar.get_nchan()
        freqs[gpu] = np.array([i0.get_Profile(0,k).get_centre_frequency() for k in range(nch)])
    ncyc = nch/npfb
    df = np.abs(np.diff(freqs[freqs.keys()[0]][:2]))[0]

    fset,edges,gpusets,gpumap,idxmap = getIndexMap(freqs, npfb, ncyc, df)

    
    return freqs,fset,gpumap,idxmap
        
def getIndexMap(freqs,npfb,ncyc,df):
    gpus = freqs.keys()
    edges = dict([(gpu,(fs[0],fs[-1])) for gpu,fs in freqs.items()])
    fset = np.array(list(set.union(*tuple([set(x) for x in freqs.values()]))))
    fset.sort()
    gpumap = np.zeros(fset.shape,dtype='int')
    idxmap = np.zeros(fset.shape,dtype='int')
    
    gpusets = []
    done = []
    for gpu in gpus:
        if gpu in done:
            continue
        done.append(gpu)
        f00,f01 = edges[gpu]
        for gpu2 in gpus:
            if gpu2 in done:
                continue
            f10,f11 = edges[gpu2]
            if np.abs(f00-f10) < ncyc*df and np.abs(f01-f11) < ncyc*df:
                gpusets.append((gpu,gpu2))
                done.append(gpu2)
                break
        else:
            gpusets.append([gpu])
            
    for gpuset in gpusets:
        if len(gpuset) == 1:
            gpu = gpuset[0]
            idx0 = np.argmin(np.abs(fset-edges[gpu][0]))
            idx1 = np.argmin(np.abs(fset-edges[gpu][1]))
            if freqs[gpu][0] < freqs[gpu][1]:
                gpumap[idx0:idx1+1] = gpu
                idxmap[idx0:idx1+1] = np.arange(idx1-idx0+1)
            else:
                gpumap[idx1:idx0+1] = gpu
                idxmap[idx1:idx0+1] = np.arange(freqs[gpu].shape[0])[::-1]
        else:
            for ipfb in range(npfb):
                for gpu in gpuset:
                    idx0 = ipfb*ncyc + ncyc/4
                    idx1 = idx0 + ncyc/2
                    f0 = freqs[gpu][idx0]
                    f1 = freqs[gpu][idx1]
                    odx0 = np.argmin(np.abs(fset-f0))
                    odx1 = np.argmin(np.abs(fset-f1))
                    print odx0,odx1
                    if odx0 < odx1:
                        gpumap[odx0:odx1] = gpu
                        idxmap[odx0:odx1] = np.arange(idx0,idx1)
                    else:
                        gpumap[odx1:odx0] = gpu
                        idxmap[odx1:odx0] = np.arange(idx0,idx1)[::-1]
    if np.any([len(x)>1 for x in gpusets]):
        # if any gpusets have more than one gpu, we have overlapping gpus, so fill in the blanks at the edges of the spectrum
        firstgpu = gpumap[ncyc/4]
        idx0 = np.argmin(np.abs(fset[0]-freqs[firstgpu]))
        idx1 = np.argmin(np.abs(fset[ncyc/4]-freqs[firstgpu]))
        gpumap[0:ncyc/4] = firstgpu
        if idx0 < idx1:
            idxmap[0:ncyc/4] = np.arange(idx0,idx1)
        else:
            idxmap[0:ncyc/4] = np.arange(idx1+1,idx0+1)[::-1]
        lastgpu = gpumap[-ncyc/4-2]
        idx0 = np.argmin(np.abs(fset[-1]-freqs[lastgpu]))
        idx1 = np.argmin(np.abs(fset[-ncyc/4-1]-freqs[lastgpu]))
        gpumap[-ncyc/4:] = lastgpu
        if idx0 < idx1:
            idxmap[-(idx1-idx0):] = np.arange(idx0,idx1)[::-1]  # not sure why this needed to be reversed...
        else:
            idxmap[-(idx0-idx1):] = np.arange(idx1+1,idx0+1)[::-1]
    return fset,edges,gpusets,gpumap,idxmap

def loadDynSpecFromCycSpecScan(scanno,mjd,gpus=None,telescope='gbt',npfb=32,thresh=0.2,alpha=-4.0,discard=0.001):
    files = getFilesForScan(scanno, mjd, telescope,gpus=gpus)
    if not files:
        raise Exception("no files found")
    fset,on,off,epochs,tints,profile,onp,offp,gpumap,idxmap,freqs,source,telescope = getOnOffSpectra(files,npfb=npfb,thresh=thresh)
    ds = dynspec.DynSpec(source=source,
                         epochs=epochs,
                         telescope=telescope,
                         tints=tints,
                         profile=profile,
                         onp=onp,
                         offp=offp,
                         freqs=fset,
                         on=on,
                         off=off,
                         guppi=False,alpha=alpha,discard=discard)
    return ds

def getOnOffSpectra(files,profile=None,npfb=32,thresh=0.2):
    #del files[7]
#    for gpu in files.keys():
#        files[gpu] = [files[gpu][0]]
#    print files
    if profile is None:
        profile = getProfile(files)
    onp = profile > np.median(profile) + (np.max(profile)-np.median(profile))*thresh
    offp = profile < np.median(profile) + (np.max(profile)-np.median(profile))*0.1
    freqs,fset,gpumap,idxmap = getFreqMap(files, npfb)
    nf = fset.shape[0]
    ons = []
    offs = []
    epochs = []
    tints = []
    for gpu in files.keys():
        osubint = 0
        oidx = np.flatnonzero(gpumap==gpu)
        iidx = idxmap[oidx]
        for fn in files[gpu]:
            print "loading: ",fn
            ar = psrchive.Archive_load(fn)
            source = ar.get_source()
            telescope = ar.get_telescope()
            
#            ar.bscrunch_to_nbin(64)
            d = ar.get_data()
#            imsk = np.zeros((d.shape[2],),dtype='bool')
#            imsk[iidx] = True
#            omsk = np.zeros((nf,),dtype='bool')
#            omsk[oidx] = True
            for k in range(d.shape[0]):
                if len(epochs) <= osubint:
                    sub = ar.get_Integration(k)
                    epochs.append(sub.get_epoch().in_days())
                    tints.append(sub.get_duration())
                try:
                    on = ons[osubint]
                    off = offs[osubint]
                except:
                    on = np.zeros((1,2,nf))
                    off = np.zeros((1,2,nf))
                    ons.append(on)
                    offs.append(off)
                x = d[k,:2,:,:]
                x = x[:,iidx,:]
                x = x[:,:,onp].mean(2)
                on[:,:,oidx] = x
                x = d[k,:2,:,:]
                x = x[:,iidx,:]
                x = x[:,:,offp].mean(2)  
                off[:,:,oidx] = x
                osubint += 1
    on = np.concatenate(ons,axis=0)
    off = np.concatenate(offs,axis=0)
    epochs = np.array(epochs)
    tints = np.array(tints)
    return fset,on,off,epochs,tints,profile,onp,offp,gpumap,idxmap,freqs,source,telescope
    
    
def getFilesForScan(scanno,mjd,telescope,gpus = None):
    filepattern = FILE_PATTERNS[telescope.lower()]
    files = {}
    if gpus is None:
        gpus = range(1,10)
    for gpu in gpus:
        g = glob.glob(filepattern % (gpu,mjd,scanno))
        if g:
            g.sort()
            files[gpu] = g
    return files

def fixCyclicDedisp(fname, nchan=32, overwrite=False, ext='fix'):
    # copied from paul's fix_cyclic_dedisp script
    import psrchive
    import os
    import psr_utils
    arch = psrchive.Archive_load(fname)
    cf = arch.get_centre_frequency()
    bw = arch.get_bandwidth()
    f_lo = cf - bw/2.0
    nch = arch.get_nchan()
    pfb_nch = nchan
    pfb_bw = bw / pfb_nch
    chan_per_pfb = nch / pfb_nch
    dm = arch.get_dispersion_measure()
    for isub in range(arch.get_nsubint()):
        sub = arch[isub]
        per = sub.get_folding_period()
        for ichan in range(nch):
            pfb_cf = f_lo + ((ichan/chan_per_pfb)+0.5)*pfb_bw
            dt = psr_utils.delay_from_DM(dm,pfb_cf) - psr_utils.delay_from_DM(dm,cf)
            for ipol in range(sub.get_npol()):
                prof = sub.get_Profile(ipol,ichan)
                prof.rotate_phase(dt/per)
    #arch.set_dedispersed(True) # doesn't work, lame...
    if (overwrite):
        outf = fname
    else:
        outf = fname + '.' + ext
    arch.unload(outf)
    os.system("psredit -m -c dmc=1 %s" % outf)
