import dynspec
import numpy as np
import glob

def getEstimates():
    pkls = glob.glob('/lakitu/scratch/nanograv/cs/analysis/dynspec/pkls/*.pkl')
    fdif = {}
    tdif = {}
    for pkl in pkls:
        dset = dynspec.unpickle(pkl)
        source = dset.source
        if fdif.has_key(source):
            fdif[source].append(dset.fit.params['fdif'].value)
            tdif[source].append(dset.fit.params['tdif'].value)
        else:
            fdif[source] = [dset.fit.params['fdif'].value]
            tdif[source] = [dset.fit.params['tdif'].value]
    return fdif,tdif
def printEstimates(fdif,tdif):
    sources = fdif.keys()
    sources.sort()
    print "fdif"
    print "="*80
    for k in sources:
        print k," : ",fdif[k],','
    
    print "tdif"
    print "="*80
        
    for k in sources:
        print k, " : ",tdif[k],','

fdifs = {
('0023+0923',430)  :  0.6,
('0030+0451',430)  :  15.0 ,
('0613-0200',800)  :  .12,
('0613-0200',1470)  :  4.0,
('1012+5307',800)  :  3,
('1455-3330',800)  :  20,
('1600-3053',1400)  :  0.2,
('1600-3053',800)  :  0.03,
('1614-2230',800)  :  0.25,
('1640+2224',430)  :  0.5,
('1643-1224',1400)  : 0.015,
('1643-1224',800)  : 0.002,
('1713+07',1400)  :  10.0,
('1713+07',800)  :  2.0,
('1741+1351',430)  :  0.2,
('1744-1134',800)  :  2.0,
('1744-1134',1400)  :  40,
('1853+1303',430)  :  0.1,
('1855+09',430)  :  0.03,
('1909-3744',800)  : 4.0,
('1918-0642',800)  :  2.0,
('1918-0642',1200)  :  8.0,
('1923+2515',430)  :  0.3,
('1944+0907',430)  :  0.06,
#2017+0603  :  [0.085788889241824817] ,  # RFI dominated
('2145-0750',1400)  :  20.0,
('2145-0750',800)  :  4.0,
('2317+1439',327)  :  0.1,
('2317+1439',430)  :  0.4,
('B1937+21', 327)  :  0.004,
('B1937+21', 430)  :  0.007,
('B1937+21', 800)  :  0.01,
('B1937+21', 1200)  :  0.25,
('B1937+21', 1400)  :  0.6,
('B1953+29',327)  :  0.001,
('B1953+29',430)  :  0.002,
('B1953+29',1460)  :  0.08,
('J0340+4130',800)  :  0.1,
('J0340+4130',1470)  :  5.0,
('J0645+51',800)  :  40.0,
('J0931-1902',800)  :  2.0,
('J0931-1902',1400)  :  40.0,
('J2010-1323',800)  :  .25,
('J2043+1711',430)  :  0.5,
('J2229+2643',430)  :  0.4,
('J2229+2643',327)  :  0.06,
('J2234+09',430)  :  1.0,
('J2302+4442',800)  :  0.2,
}
tdifs= {
('0023+0923',430)  :  1000,
('0030+0451',430)  :  2000 ,
('0613-0200',800)  :  350,
('0613-0200',1470)  :  600,
('1455-3330',800)  :  1000,
('1600-3053',1400)  :  400,
('1600-3053',800)  :  200,
('1614-2230',800)  :  400,
('1640+2224',430)  :  1000,
('1643-1224',1400)  : 70,
('1643-1224',800)  : 30,
('1713+07',1400)  :  1600.0,
('1713+07',800)  :  600.0,
('1741+1351',430)  :  700,
('1744-1134',800)  :  700.0,
('1744-1134',1400)  :  1500,
('1853+1303',430)  :  800,
('1855+09',430)  :  450,
('1909-3744',800)  : 1500.0,
('1918-0642',800)  :  400.0,
('1918-0642',1200)  :  800.0,
('1923+2515',430)  :  650,
('1944+0907',430)  :  150,
#2017+0603  :  [0.085788889241824817] ,  # RFI dominated
('2145-0750',1400)  :  1500.0,
('2145-0750',800)  :  600.0,
('2317+1439',327)  :  250,
('2317+1439',430)  :  750,
('B1937+21', 327)  :  70,
('B1937+21', 430)  :  70,
('B1937+21', 800)  :  150,
('B1937+21', 1200)  :  150,
('B1937+21', 1400)  :  150,
('B1953+29',327)  :  50,
('B1953+29',430)  :  50,
('B1953+29',1460)  :  100,
('J0340+4130',800)  :  200,
('J0340+4130',1400)  :  600,
#('J0931-1902',1400)  :  10.0,
('J2010-1323',800)  :  600,
('J2043+1711',430)  :  1000,
('J2229+2643',430)  :  500,
('J2229+2643',327)  :  400,
('J2234+09',430)  :  500,
('J2302+4442',800)  :  600,
}

fdif_by_band = {}
tdif_by_band = {}
for ((source,band),val) in fdifs.items():
    if fdif_by_band.has_key(source):
        fdif_by_band[source][band] = val
    else:
        fdif_by_band[source] = {band:val}

for ((source,band),val) in tdifs.items():
    if tdif_by_band.has_key(source):
        tdif_by_band[source][band] = val
    else:
        tdif_by_band[source] = {band:val}
        

def _findEstimate(source,freq,bw,index):
    if index.has_key(source):
        bands = np.array(index[source].keys())
        band = bands[np.abs(bands-freq) < bw]
        if len(band):
            return index[source][band[0]]
    return None

def findEstimates(source,freq,bw):
    tdif0 = _findEstimate(source, freq, bw, tdif_by_band)
    fdif0 = _findEstimate(source, freq, bw, fdif_by_band)
    return tdif0,fdif0
    