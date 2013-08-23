import matplotlib
matplotlib.use('agg')
import uppispec
import dynspec
import multiprocessing
import functools
import glob
import os
from matplotlib.figure import Figure                         
from matplotlib.backends.backend_agg import FigureCanvasAgg 

telescope='gbt'
gpus=None#[1,2,3,4,6]
redo=False
outdir = '/lakitu/data/scratch/nanograv/uppi/analysis/dynspec'
plotdir= '/home/gjones/uppi/analysis/dynspec/plots'
#todo = glob.glob('/lakitu/data/scratch/nanograv/puppi/*/*uppi_?????_*_*_0001.fits')
todo = glob.glob('/lakitu/data/scratch/nanograv/guppi/*/*uppi_?????_*_*_0001.fits')
todo = [x for x in todo if x.find('cal') < 0]
todo.sort()
#todo = todo[:2]
done = set(glob.glob(os.path.join(outdir,('pkls/ds_*_%s_*.pkl' % telescope))))
errors = {}
def pOne(fn):
    return uppispec.processOne(fn,done,telescope,outdir,plotdir,redo)
p = multiprocessing.Pool(12)
errors = dict(p.map(pOne,todo))
dynspec.pickle('uppierrors.pkl', errors)
print errors