
import numpy
import datetime
import matplotlib.pyplot
import astropy.io.fits
import astropy.coordinates
from astropy.units import deg, arcmin, arcsec, MHz, Hz

import bcat.common.common_structure


def create_stage1(record):
    pass
    

def get_obstime(record):
    def decode(timestamp):
        t = datetime.datetime.strptime(timestamp, '%Y-%m-%dT%H:%M:%S').timestamp()
        return t - 9

    obstime = astropy.time.Time(numpy.array([decode(_t) for _t in record['date-obs']]), format='unix')
    return obstime


def get_rf_header(record):
    lo1 = record['lofreq'] * MHz
    lo2 = record['_2ndlo'] * Hz
    sideband1 = record['sideband']
    sideband2 = record['_2ndsb']

    sb1 = numpy.zeros_like(sideband1, dtype='int8')
    sb1[sideband1 == 'U'] = 1
    sb1[sideband1 == 'L'] = -1

    sb2 = numpy.zeros_like(sideband2, dtype='int8')
    sb2[sideband2 == 'U'] = 1
    sb2[sideband2 == 'L'] = -1

    crval = (lo1 + (sb1 * lo2)).to('Hz').value
    crpix = numpy.zeros_like(crval)
    cdelt = sb2 * record['freqres']
    
    return bcat.common.common_structure.freq_axis_array(cdelt, crpix, cdelt)
    
    
def get_coord(record):
    x = record['crval2'] * deg + record['lamdel'] * arcsec
    y = record['crval3'] * deg + record['betdel'] * arcsec

    if record['coordsys'][0].lower() == 'b1950':
        frame = 'fk4'
    elif record['coordsys'][0].lower() == 'j2000':
        frame = 'fk5'
    elif record['coordsys'][0].lower() == 'galactic':
        frame = 'galactic'
        pass

    obstime = get_obstime(record)
    coord = astropy.coordinates.SkyCoord(x, y, frame=frame, obstime=obstime)
    return coord

    

def print_header(record):
    def register(key, description):
        item = record[key.lower()][0]
        
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
        register('OBJECT', 'Obs object name'),
        register('BANDWID', 'Bandwidth of spectrum'),
        register('DATE-OBS', 'Obs start timestamp'),
        register('TDIM6', 'Shape of each record'),
        register('TUNIT6', 'Unit of spectrum'),
        register('CTYPE1', 'z : type'),
        register('CRVAL1', 'z : ref val'),
        register('CRPIX1', 'z : ref pix'),
        register('CDELT1', 'z : delta'),
        register('CTYPE2', 'x : type'),
        register('CRVAL2', 'x : ref val'),
        register('CTYPE3', 'y : type'),
        register('CRVAL3', 'y : ref val'),
        register('OBSERVER', 'Observer name'),
        register('OBSMODE', 'Observation mode'),
        register('MOLECULE', 'Target molecule'),
        register('TRANSITI', 'Target transition'),
        register('FRONTEND', 'Receiver name'),
        register('BACKEND', 'Backend name'),
        register('FREQRES', 'Frequency resolution'),
        register('TIMESYS', ''),
        register('VELDEF', 'Definition of velocity'),
        register('RESTFREQ', 'Rest freq. of line'),
        register('COORDSYS', ''),
        register('COSYDEL', ''),
        register('LOFREQ', ''),
        register('SYNTH', ''),
        register('SIDEBAND', ''),
        register('_2NDSB', ''),
        register('_3RDSB', ''),
        register('_2NDLO', ''),
        register('_3RDLO', ''),
    ]

    for _ in txt:
        print(f'{_[1]:22s} : {_[0]:8s} : {_[2]}')
        continue

    return 


def plot_header(record, title=''):
    def _plot(ax, key):
        y = record[key]
        ax.plot(y)
        ax.set_xticklabels('')
        ax.set_title(key.upper())
        return
    
    keys = [
        'lamdel',
        'betdel',
        'azimuth',
        'elevatio',
        'crval2',
        'crval3',
        'exposure',
        'thot',
        'tcold',
        'tambient',
        'pressure',
        'humidity',
        'otfvlam',
        'otfvbet',
        'otfscann',
        'subscan',
        'vframe',
        'vframe2',
        'synth',
        '_2ndlo',
    ]
    
    fig = matplotlib.pyplot.figure(figsize=(18, 8))
    ax = [fig.add_subplot(4, 6, i) for i in range(1, 21)]
    [_plot(_ax, _key) for _ax, _key in zip(ax, keys)]
    fig.suptitle(f'{title}')
    return fig

def plot_spectra(record):
    sp = record['data']
    c = numpy.median(sp)
    std = numpy.std(sp)
    vmin = c - 1 * std
    vmax = c + 2 * std
    
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111, aspect='auto')
    ax.imshow(sp, vmin=vmin, vmax=vmax, origin='lower', aspect='auto')
    ax.set_xlabel('velocity (ch)')
    ax.set_ylabel('# of spectra')
    ax.set_title(f'')
    return fig
