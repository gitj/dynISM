import matplotlib
matplotlib.use('agg')
import cycspec
import dynspec
import multiprocessing
import functools
import glob
import os
from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 

outdir = '/lakitu/data/scratch/nanograv/cs/analysis/dynspec'
plotdir= '/home/gjones/cs/analysis/dynspec/rfiplots'
todo = glob.glob(os.path.join(outdir,'pkls/ds*.pkl'))
todo.sort()
errors = {}
repickle=True
def pOne(fn):
    try:
        fig = Figure(figsize=(10,12))
        fig.subplots_adjust(left=0.09,bottom=0.05,top=0.95,right=0.95)
        ds = dynspec.unpickle(fn)
        ds = dynspec.DynSpec(ds.source,ds.epochs,ds.telescope,ds.tints,ds.profile,ds.onp,ds.offp,ds.freqs,ds.on,ds.off,ds.guppi,ds.alpha,discard=0.05)
        ds.refit()
        ds.plot(fig=fig)
        outfn = os.path.split(fn)[1]
        plotname = os.path.join(plotdir,outfn + '.png')
        esc_fname = outfn.replace('_',r'\_')
        fig.suptitle(('%s @ %s %s' % (ds.source,ds.telescope,esc_fname)),size='medium')
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(plotname)
        if repickle:
            dynspec.pickle(fn, ds)
    except Exception, e:
        return fn,e
    return fn,None
p = multiprocessing.Pool(12)
errors = dict(p.map(pOne,todo))
cycspec.dynspec.pickle('plot_errors.pkl', errors)
print errors