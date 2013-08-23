import matplotlib
matplotlib.use('agg')
import cycspec
import multiprocessing
import functools
import glob
import os

telescope='ao'
gpus=None#[1,2,3,4,6]
redo=False
outdir = '/lakitu/data/scratch/nanograv/cs/analysis/dynspec'
plotdir= '/home/gjones/cs/analysis/dynspec/plots'
todo = glob.glob((os.path.split(cycspec.FILE_PATTERNS[telescope])[0] +'/*_nth4_0001.ar.fix')% 1)
todo.sort()
#todo = todo[:2]
done = set(glob.glob(os.path.join(outdir,('pkls/ds_*_%s.pkl' % telescope))))
errors = {}
def pOne(fn):
    return cycspec.processOne(fn, done, gpus, telescope, outdir, plotdir, redo)
p = multiprocessing.Pool(12)
errors = dict(p.map(pOne,todo))
cycspec.dynspec.pickle('errors.pkl', errors)
print errors