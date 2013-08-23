import lmfit
import numpy as np
import functools

def gammaIs(params,x):
    a = params['tdif'].value
    return np.exp(-(x/a)**(5/3.0))

def gammaIs2(params,x):
    a = params['tdif'].value
    offs = params['offs'].value
    return (np.exp(-(x/a)**(5/3.0))+offs)/(offs+1.0)

def gammaIs3(params,x):
    a = params['tdif'].value
    offs = params['offs'].value
    scale = params['scale'].value
    return scale*np.exp(-(x/a)**(5/3.0))+offs

def gammaIf(params,x):        
    a = params['fdif'].value
    return np.exp(-((np.log(2)*x/a)))

def gammaIf2(params,x):        
    a = params['fdif'].value
    scale = params['scale'].value
    return scale*np.exp(-((np.log(2)*x/a)))

def gammaIf3(params,x):        
    a = params['fdif'].value
    scale = params['scale'].value
    offs = params['offs'].value
    return scale*np.exp(-((np.log(2)*x/a)))+offs

def fitIs(x,y):
    params = lmfit.Parameters()
    params.add('tdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()])
    return Fit(x,y,gammaIs,params)

def fitIs2(x,y):
    params = lmfit.Parameters()
    params.add('tdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()])
    params.add('offs',value=0.5)
    return Fit(x,y,gammaIs2,params)

def fitIs3(x,y):
    params = lmfit.Parameters()
    params.add('tdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()])
    params.add('offs',value=y.min())
    params.add('scale',value=y.max())
    return Fit(x,y,gammaIs3,params)

def fitIf(x,y):
    params = lmfit.Parameters()
    params.add('fdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()], min = 0.0)
    return Fit(x,y,gammaIf,params)

def fitIf2(x,y):
    params = lmfit.Parameters()
    params.add('fdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()], min = 0.0)
    params.add('scale',value=y.max(),min=0.0)
    return Fit(x,y,gammaIf2,params)

def fitIf3(x,y,fdif0=None):
    params = lmfit.Parameters()
    if fdif0 is None:
        params.add('fdif',value = x[np.abs(y-(y.max()-y.min())/2.0).argmin()],min=0,max=400.0)
    else:
        params.add('fdif',value = fdif0,min=fdif0/5.0,max=fdif0*5.0)
    params.add('offs',value=y.min())
    params.add('scale',value=y.max())
    return Fit(x,y,gammaIf3,params)


def simFitAcf(t,tacf,f,facf,tdif0=None,fdif0=None):
    params = lmfit.Parameters()
    if tdif0 is None:
        params.add('tdif',value = t[np.abs(tacf-(tacf.max()-tacf.min())/2.0).argmin()], min = 0.0,max=20000.0)
    else:
        params.add('tdif',value = tdif0, min = tdif0/5.0,max=tdif0*5.0)
    if fdif0 is None:
        params.add('fdif',value = f[np.abs(facf-(facf.max()-facf.min())/2.0).argmin()], min = 0.0,max=400.0)
    else:
        params.add('fdif',value = fdif0, min = fdif0/5.0, max=fdif0*5.0)
    params.add('scale',value=tacf.max(),min=0.0)
    params.add('offs',value=facf.min(),min=0.0)
    def resid(p,tt,ttacf,ff,ffacf):
        fresid = ffacf - gammaIf3(p, ff)
        tresid = ttacf - gammaIs3(p, tt)
        return np.hstack((fresid,tresid))
    
    mi = lmfit.minimize(resid,params,args=(t,tacf,f,facf))
    return mi

def makeResid(func):
    return functools.partial(genericResid,func)
def genericResid(func,params,x,y):
    return y - func(params,x)



class Fit():
    def __init__(self,x,y,func,params,args=()):
        self.params = params
        self.func = func
        self.rfunc = makeResid(func)
        
        self.minimizer = lmfit.minimize(self.rfunc, params, args=(x,y)+args)
        
    