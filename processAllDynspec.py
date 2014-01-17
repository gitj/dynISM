import matplotlib
matplotlib.use('agg')
import cycspec
import multiprocessing
import functools
import glob
import os
from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 
import dynspec

telescope='ao'
gpus=None#[1,2,3,4,6]
redo=True #False
outdir = '/lakitu/data/scratch/nanograv/cs/analysis/dynspec'
plotdir= '/home/gjones/cs/analysis/dynspec/plots'
rfiplotdir= '/home/gjones/cs/analysis/dynspec/rfiplots'
todo = glob.glob((os.path.split(cycspec.FILE_PATTERNS[telescope])[0] +'/*_nth4_0001.ar.fix')% 1)
todo.sort()
#todo = todo[:2]
done = set(glob.glob(os.path.join(outdir,('pkls/ds_*_%s.pkl' % telescope))))
errors = {}
repickle=True
def pOne(fn):
    fn, e,pkl = cycspec.processOne(fn, done, gpus, telescope, outdir, plotdir, redo)
    if e:
        return fn, e
    pkl = os.path.join(outdir,'pkls',pkl)
    print pkl
    try:
        fig = Figure(figsize=(10,12))
        fig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
        ds = dynspec.unpickle(pkl)
        ds = dynspec.DynSpec(ds.source,ds.epochs,ds.telescope,ds.tints,ds.profile,ds.onp,ds.offp,ds.freqs,ds.on,ds.off,ds.guppi,ds.alpha,discard=0.05)
        ds.refit()
        ds.plot(fig=fig)
        outfn = os.path.split(pkl)[1]
        plotname = os.path.join(rfiplotdir,outfn + '.png')
        esc_fname = outfn.replace('_',r'\_')
        fig.suptitle(('%s @ %s %s' % (ds.source,ds.telescope,esc_fname)),size='medium')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(plotname)
        if repickle:
            dynspec.pickle(pkl, ds)
    except Exception, e:
        return fn,e
    return fn,None

#def pOne(fn):
#    return cycspec.processOne(fn, done, gpus, telescope, outdir, plotdir, redo)
p = multiprocessing.Pool(1)
errors = dict(p.map(pOne,todo))
cycspec.dynspec.pickle('errors.pkl', errors)
print errors