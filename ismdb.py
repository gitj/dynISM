import gluon
from gluon.dal import DAL,Field

ismdb = DAL('sqlite://nano.db', folder = '/home/gjones/workspace/dynism/tables')

Instruments = ismdb.define_table('instruments',
                              Field('name'),
                              Field('telescope'))
Sources = ismdb.define_table('sources',
                          Field('name'))
Observations = ismdb.define_table('observations',
                               Field('mjd','integer'),
                               Field('scan','integer'),
                               Field('source',ismdb.sources),
                               Field('epoch','double'),
                               Field('ref_freq','double'),
                               Field('instrument',ismdb.instruments),
                               Field('total_duration','double'),
                               Field('subint_duration','double'),
                               Field('bandwidth','double'),
                               Field('nchan','integer'),
                               Field('tdif_acf','double'),
                               Field('tdif_acf_err','double'),
                               Field('fdif_acf','double'),
                               Field('fdif_acf_err','double'),
                               Field('alpha','double'),
                               Field('plot_file','string')
                               )
Flags = ismdb.define_table('flags',
                        Field('mjd','integer'),
                        Field('scan','integer'),
                        Field('instrument',ismdb.instruments),
                        Field('flag_type',requires=gluon.IS_IN_SET(['tdif','fdif']),
                              readable=False,writable=False),
                        Field('flag_value',label='Flag'),
                        )

def get_band(freq):
    if 270 < freq < 370:
        return '327'
    if 380 < freq < 450:
        return '430' 
    if 700 < freq < 950:
        return '820'
    if 1000 < freq < 1800:
        return 'Lband'
    if 1800 < freq < 2500:
        return 'Sband'
    return 'unknown %.1f' % freq

band_edges = {'327': (270,370),
              '430': (380,450),
              '820': (700,950),
              'Lband': (1000,1800),
              'Sband': (1800,2500)}
def normalize_source(raw):
    aliases = {'1713+07':'1713+0747',
               '0645+5150':'0645+5158',
               '0645+51':'0645+5158',
               '1640+22':'1640+2224',
               '1923+25' : '1923+2515',
               '1955+2908':'1953+29'
               }
    if raw[0].isalpha():
        raw = raw[1:]
    if raw in aliases:
        raw = aliases[raw]
    return raw

def populate():
    import dynspec
    import glob
    import re
    import numpy as np
    import os
    fndecode = re.compile(r".*ds_(?P<source>[BJ]*\d{4}[+-]\d{2,4})_(?P<mjd>\d*)_(?P<scan>\d*)_(?P<telescope>.*)_(?P<instr>.*)\.pkl")
    pklglob = glob.glob('/lakitu/scratch/nanograv/uppi/analysis/dynspec/pkls/*.pkl')
    plot_dir = '/home/gjones/uppi/analysis/dynspec/plots'
    
    pklglob = glob.glob('/lakitu/scratch/nanograv/cs/analysis/dynspec/pkls/ds*gbt.pkl')+ glob.glob('/lakitu/scratch/nanograv/cs/analysis/dynspec/pkls/ds*ao.pkl')
    plot_dir = '/home/gjones/cs/analysis/dynspec/rfiplots'
    fndecode = re.compile(r".*ds_(?P<source>[BJ]*\d{4}[+-]\d{2,4})_(?P<mjd>\d*)_(?P<scan>\d*)_(?P<telescope>.*)\.pkl")
    for pkl in pklglob:
        print pkl
        info = fndecode.search(pkl).groupdict()
        ds = dynspec.unpickle(pkl)
        normsource = normalize_source(ds.source)
        source_id = ismdb(ismdb.sources.name == normsource).select('id')
        if not source_id:
            print "inserting source"
            source_id = ismdb.sources.insert(name=normsource)
            ismdb.commit()
        else:
            source_id = source_id[0]
        inst_map = {(True,'gbt') : 'guppi',
                    (True,'arecibo') : 'puppi',
                    (False,'gbt') : 'cycspec-g',
                    (False,'arecibo') : 'cycspec-a'}
        inst_name = inst_map[(ds.guppi,ds.telescope.lower())]
        inst_id = ismdb(ismdb.instruments.name == inst_name).select('id')
        if inst_id:
            inst_id = inst_id[0]
        else:
            print "inserting instrument"
            inst_id = ismdb.instruments.insert(name=inst_name,telescope=ds.telescope.lower())
            ismdb.commit()
        if ds.guppi:
            ref_freq = ds.fc
            tdif = ds.fit_stretch.params['tdif'].value
            tdif_err = ds.fit_stretch.params['tdif'].stderr
            fdif = ds.fit_stretch.params['fdif'].value
            fdif_err = ds.fit_stretch.params['fdif'].stderr
            plot_file = os.path.join(plot_dir,('ds_%s_%s_%s_%s_%s.pkl.png' % 
                                               (info['source'],info['mjd'],info['scan'],
                                                info['telescope'],info['instr'])))
        else:
            ref_freq = ds.freqs.mean()
            tdif = ds.fit.params['tdif'].value
            tdif_err = ds.fit.params['tdif'].stderr
            fdif = ds.fit.params['fdif'].value
            fdif_err = ds.fit.params['fdif'].stderr
            plot_file = os.path.join(plot_dir,('ds_%s_%s_%s_%s.pkl.png' % 
                                               (info['source'],info['mjd'],info['scan'],
                                                info['telescope'])))
        ismdb.observations.update_or_insert((ismdb.observations.mjd == int(info['mjd'])) & (ismdb.observations.scan == int(info['scan'])),
                               mjd = int(info['mjd']),
                               scan = int(info['scan']),
                               source = source_id,
                               instrument = inst_id,
                               epoch = ds.epochs.mean(),
                               ref_freq = ref_freq,
                               total_duration = ds.tints.sum(),
                               subint_duration = np.median(ds.tints),
                               bandwidth = ds.freqs.ptp(),
                               nchan = ds.freqs.shape[0],
                               alpha = ds.alpha,
                               tdif_acf = tdif,
                               tdif_acf_err = tdif_err,
                               fdif_acf = fdif,
                               fdif_acf_err = fdif_err,
                               plot_file = plot_file
                               )
        ismdb.commit()