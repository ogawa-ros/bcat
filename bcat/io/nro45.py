
import io
import struct
import datetime
import pathlib

import numpy
import matplotlib.pyplot
import astropy.time
import astropy.coordinates
import astropy.units

import bcat


size_obs_header = 15136
size_scan_record = 6568

dtype_obs_header = [
    ('lofil0', 'S8'),
    ('ver0', 'S8'),
    ('group0', 'S16'),
    ('proj0', 'S16'),
    ('sched0', 'S24'),
    ('obsvr0', 'S40'),
    ('lostm0', 'S16'),
    ('loetm0', 'S16'),
    ('arynm0', '>i4'),
    ('nscan0', '>i4'),
    ('title0', 'S120'),
    ('obj0', 'S16'),
    ('epoch0', 'S8'),
    ('ra00', '>f8'),
    ('dec00', '>f8'),
    ('gl00', '>f8'),
    ('gb00', '>f8'),
    ('ncalb0', '>i4'),
    ('scncd0', '>i4'),
    ('scmod0', 'S120'),
    ('vel0', '>f8'),
    ('vref0', 'S4'),
    ('vdef0', 'S4'),
    ('swmod0', 'S8'),
    ('frqsw0', '>f8'),
    ('dbeam0', '>f8'),
    ('mltof0', '>f8'),
    ('cmtq0', '>f8'),
    ('cmte0', '>f8'),
    ('cmtsom0', '>f8'),
    ('cmtnode0', '>f8'),
    ('cmti0', '>f8'),
    ('cmttm0', 'S24'),
    ('sbdx0', '>f8'),
    ('sbdy0', '>f8'),
    ('sbdz10', '>f8'),
    ('sbdz20', '>f8'),
    ('dazp0', '>f8'),
    ('delp0', '>f8'),
    ('cbind0', '>i4'),
    ('nch0', '>i4'),
    ('chrange0', '>2i4'),
    ('alctm0', '>f8'),
    ('iptim0', '>f8'),
    ('pa0', '>f8'),
    ('rx0', '35S16'),
    ('hpbw0', '>35f8'),
    ('effa0', '>35f8'),
    ('effb0', '>35f8'),
    ('effl0', '>35f8'),
    ('efss0', '>35f8'),
    ('gain0', '>35f8'),
    ('horn0', '35S4'),
    ('poltp0', '35S4'),
    ('poldr0', '>35f8'),
    ('polan0', '>35f8'),
    ('dfrq0', '>35f8'),
    ('sidbd0', '35S4'),
    ('refn0', '>35i4'),
    ('ipint0', '>35i4'),
    ('multn0', '>35i4'),
    ('mltscf0', '>35f8'),
    ('lagwin0', '35S8'),
    ('bebw0', '>35f8'),
    ('beres0', '>35f8'),
    ('chwid0', '>35f8'),
    ('arry0', '>35i4'),
    ('nfcal0', '>35i4'),
    ('f0cal0', '>35f8'),
    ('fqcal0', '>(35,10)f8'),
    ('chcal0', '>(35,10)f8'),
    ('cwcal0', '>(35,10)f8'),
    ('scnlen0', '>i4'),
    ('sbind0', '>i4'),
    ('bit0', '>i4'),
    ('cdmy1', 'V188'),
]

dtype_scan_record = [
    ('lsfil0', 'S4'),
    ('iscn0', '>i4'),
    ('lavst0', 'S24'),
    ('scntp0', 'S8'),
    ('dscx0', '>f8'),
    ('dscy0', '>f8'),
    ('scx0', '>f8'),
    ('scy0', '>f8'),
    ('paz0', '>f8'),
    ('pel0', '>f8'),
    ('raz0', '>f8'),
    ('rel0', '>f8'),
    ('xx0', '>f8'),
    ('yy0', '>f8'),
    ('arryt0', 'S4'),
    ('temp0', '>f4'),
    ('patm0', '>f4'),
    ('ph200', '>f4'),
    ('vwind0', '>f4'),
    ('dwind0', '>f4'),
    ('tau0', '>f4'),
    ('tsys0', '>f4'),
    ('batm0', '>f4'),
    ('line0', '>i4'),
    ('flag0', '>i4'),
    ('rms0', '>f8'),
    ('unknown', '4c'),
    ('vrad0', '>f8'),
    ('frq00', '>f8'),
    ('fqtrk0', '>f8'),
    ('fqif10', '>f8'),
    ('alcv0', '>f8'),
    ('offcd0', '>(2,2)f8'),
    ('idmy0', '>i4'),
    ('cdmy1', '156c'),
    ('sfctr0', '>f8'),
    ('adoff0', '>f8'),
    ('ldata', '>6144B'),
]


def get_array_name(record):
    arr = record[0]['arryt0'].decode()
    if arr in ['A1', 'A2', 'A3', 'A4', 'A5', 
               'A6', 'A7', 'A8', 'A9']:
        arr = arr[0] + '0' + arr[1]
        pass
    return arr

def convert_array_name(name):
    if name in ['A1', 'A2', 'A3', 'A4', 'A5', 
                'A6', 'A7', 'A8', 'A9']:
        name = name[0] + '0' + name[1]
        pass
    return name
    

def read_header_item(stream, data_format):
    size = struct.calcsize(data_format)
    item = struct.unpack(data_format, stream.read(size))
    
    if data_format.endswith('s'):
        item = tuple(_.decode().strip('\x00') for _ in item)
    
    elif data_format.endswith('c'):
        item = b''.join(item)
        pass
                
    if len(item) == 1:
        item = item[0]
        pass

    return item
    

def read_obs_header(filepath):
    f = open(filepath, 'rb')
    d = f.read(size_obs_header)
    return numpy.frombuffer(d, dtype=dtype_obs_header)

def read_scan_record(filepath):
    f = open(filepath, 'rb')
    f.seek(size_obs_header, 0)
    d = f.read()
    return numpy.frombuffer(d, dtype=dtype_scan_record)

def split(scan_records):
    arr = set(scan_records['arryt0'])
    return {convert_array_name(_a.decode()): scan_records[scan_records['arryt0'] == _a]
            for _a in arr}

def make_spectra(scan_records):
    def decode(raw_spec):
        raw_spec = raw_spec.astype('uint32')
        i0 = raw_spec[0::3]
        i1 = raw_spec[1::3]
        i2 = raw_spec[2::3]
        ii0 = (i0 << 4) + (i1 >> 4)
        ii1 = ((i1 & 0xF) << 8) + (i2)
        return numpy.array([ii0, ii1]).T.flatten()
    
    sp = numpy.array([decode(_r['ldata']) for _r in scan_records])
    return (sp.T * scan_records['sfctr0'] + scan_records['adoff0']).T.astype('float32')
    

def create_stage2(obs_header, scan_records, spectrum, label):
    x = scan_records['scx0'] * astropy.units.rad
    y = scan_records['scy0'] * astropy.units.rad
    
    if obs_header['scncd0'] == 0:
        if obs_header['epoch0'] == b'J2000':
            frame = 'fk5'
        elif obs_header['epoch0'] == b'B1950':
            frame = 'fk4'
            pass

    elif obs_header['scncd0'] == 1:
        frame = 'galactic'
        pass
    
    def get_obstime(timestamp):
        t = timestamp.decode()
        tt = datetime.datetime.strptime(t, '%Y%m%d%H%M%S.%f').timestamp()
        return tt
    
    obstime = astropy.time.Time(numpy.array([get_obstime(_t) for _t in scan_records['lavst0']]) - 9, format='unix')
    rf = bcat.structure.freq_axis(0 * astropy.units.Hz, 0 * astropy.units.Hz)
    coord = astropy.coordinates.SkyCoord(x, y, frame=frame, obstime=obstime)
    d = bcat.structure.stage2_data(label=label, coord=coord, spectrum=spectrum, rf=rf)
    c = bcat.stage2.container(d)
    return c


def load(filepath):
    raw = nro45_raw(filepath)
    stage2 = {key: create_stage2(raw.obs_header,
                                 raw.records[key],
                                 raw.spectra[key],
                                 f'{filepath}.{key}')
              for key in raw.records}
    return stage2

def load_as_rawdata(filepath):
    return nro45_raw(filepath)

def print_obs_header(header):
    def register(key, description):
        item = header[key.lower()][0]
        
        if isinstance(item, bytes):
            item = item.decode()
            pass

        if isinstance(item, numpy.ndarray):
            if isinstance(item[0], bytes):
                item = [_.decode() for _ in item]
                pass
            pass
            
        
        return [key, description, item]
        
    txt = [
        register('LOFIL0', 'File type'),
        register('VER0', 'File version'),
        register('GROUP0', 'Group name'),
        register('PROJ0', 'Project name'),
        register('SCHED0', 'Obs file name'),
        register('OBSVR0', 'Observer'),
        register('LOSTM0', 'Obs start timestamp'),
        register('LOETM0', 'Obs end timestamp'),
        register('ARYNM0', 'Number of array'),
        register('NSCAN0', 'Number of records'),
        register('SCNLEN0', 'Size of scan record'),
        register('TITLE0', 'Obs title'),
        register('OBJ0', 'Obs object name'),
        register('EPOCH0', 'Epoch'),
        register('RA00', 'Center R. A. (rad)'),
        register('DEC00', 'Center Dec. (rad)'),
        register('GL00', 'Center GLON (rad)'),
        register('GB00', 'Center GLAT (rad)'),
        register('NCALB0', 'Calibration interval'),
        register('SCNCD0', 'Coord (0 RADec, 1 LB, 2 AzEl)'),
        register('SCMOD0', 'Sequence pattern'),
        register('VEL0', 'Vobj (m/s)'),
        register('VREF0', 'Frame of ref. for V'),
        register('VDEF0', 'Definition for V'),
        register('SWMOD0', 'SW Mode'),
        register('FRQSW0', 'Freq. of beam SW'),
        register('DBEAM0', 'Dist. of beam SW (rad)'),
        register('MLTOF0', 'Rx rot angle (rad)'),
        register('SBDX0', 'Sub Ref dx (mm)'),
        register('SBDY0', 'Sub Ref dy (mm)'),
        register('SBDZ10', 'Sub Ref dz1 (mm)'),
        register('SBDZ10', 'Sub Ref dz2 (mm)'),
        register('DAZP0', 'Pointing daz (rad)'),
        register('DELP0', 'Pointing del (rad)'),
        register('CBIND0', 'Number of binding'),
        register('NCH0', 'Number of channel'),
        register('CHRANGE0', 'Channel range'),
        register('ALCTM0', 'ALC time constant'),
        register('IPTIM0', 'Integration time (s)'),
        register('PA0', 'Position angle (rad)'),
        register('RX0', 'Rx frontend'),
        register('SIDBD0', 'Sideband'),
        register('REFN0', 'SG number'),
        register('MULTN0', 'Beam number'),
        register('BEBW0', 'BE bandwidth (Hz)'),
        register('BERES0', 'BE resolution (Hz)'),
        register('CHWID0', 'BE ch separation'),
        register('ARRY0', 'BE array status of use'),
        register('NFCAL0', 'Num of freq calb'),
        register('F0CAL0', 'Center freq of array'),
        register('FQCAL0', ''),
        register('CHCAL0', ''),
        register('CWCAL0', ''),
    ]

    for _ in txt:
        print(f'{_[1]:22s} : {_[0]:7s} : {_[2]}')
        continue

    return 


def plot_scan_record_header(record, title=''):
    def _plot(ax, key):
        x = record['iscn0'][1:]
        y = record[key][1:]
        ax.plot(x, y)
        ax.set_xticklabels('')
        ax.set_title(key.upper())
        return
            
    keys = [
        'dscx0',
        'dscy0',
        'scx0',
        'scy0',
        'paz0',
        'pel0',
        'raz0',
        'rel0',
        'xx0',
        'yy0',
        'temp0',
        'patm0',
        'ph200',
        'vwind0',
        'dwind0',
        'tsys0',
        'vrad0',
        'frq00',
        'fqtrk0',
        'fqif10',
        'sfctr0',
        'adoff0',
    ]

    label = record[0]["lavst0"].decode()
    arr = get_array_name(record)
    
    fig = matplotlib.pyplot.figure(figsize=(18, 8))
    ax = [fig.add_subplot(4, 6, i) for i in range(1, 23)]
    [_plot(_ax, _key) for _ax, _key in zip(ax, keys)]
    fig.suptitle(f'{title} : {arr}')
    return fig

def plot_scan_record_spectra(record, sp=None):
    if sp is None:
        sp = make_spectra(record)
        pass

    sp = sp[1:]
    c = numpy.median(sp)
    std = numpy.std(sp)
    vmin = c - 1 * std
    vmax = c + 5 * std
    
    arr = get_array_name(record)
    
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.imshow(sp, vmin=vmin, vmax=vmax, origin='lower')
    ax.set_xlabel('velocity (ch)')
    ax.set_ylabel('# of spectra')
    ax.set_title(f'{record[0]["lavst0"].decode()} : {arr}')
    return fig



class nro45_raw:
    def __init__(self, path):
        self.path = pathlib.Path(path)
        self.filename = self.path.name
        self.obs_header = read_obs_header(self.path)
        all_records = read_scan_record(self.path)
        on_records = all_records[all_records['scntp0'] != b'ZERO']
        self.records = split(on_records)
        self.spectra = {k: make_spectra(self.records[k]) for k in self.records}
        pass

    def print_obs_header(self):
        return print_obs_header(self.obs_header)
    
    def plot_scan_header(self, arr='A01'):
        return plot_scan_record_header(self.records[arr], self.filename)

    def plot_spectra(self):
        if len(self.records) <= 16:
            fig = matplotlib.pyplot.figure(figsize=(12, 8))
            ax = [fig.add_subplot(4, 4, i) for i in range(1, 17)]
        else:
            fig = matplotlib.pyplot.figure(figsize=(12, 16))
            ax = [fig.add_subplot(8, 4, i) for i in range(1, 33)]
        for _a, _k in zip(ax, sorted(self.records)):
            sp = self.spectra[_k][1:]
            c = numpy.median(sp)
            std = numpy.std(sp)
            vmin = c - 1 * std
            vmax = c + 5 * std
            _a.imshow(sp, vmin=vmin, vmax=vmax, origin='lower')
            _a.set_title(_k)
            continue
        fig.suptitle(f'\n{self.filename}')
        return fig
