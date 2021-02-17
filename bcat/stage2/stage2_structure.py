import numpy
import astropy.units

import bcat.common.bcatio


dtype_base = [
    ('timestamp', 'float64'),
    ('x', 'float32'),
    ('y', 'float32'),
    ('f0', 'float32'),    
    ('df', 'float32'),    
    ('spec', '%dfloat32'),
    ('vobs', 'float32'),
    ('rms', 'float32'),
]

def generate_stage2_array(t, x, y, f0, df, spec, vobs=0, rms=0):
    num_spec, num_ch = spec.shape
    
    dtype = dtype_base.copy()
    dtype[5] = ('spec', f'{dtype[5][1]}'%(num_ch))
    
    d = numpy.array(numpy.zeros(num_spec), dtype=dtype)
    d['timestamp'] = t
    d['x'] = x
    d['y'] = y
    d['f0'] = f0
    d['df'] = df
    d['spec'] = spec
    d['rms'] = rms
    return d


def save_stage2(path, t, x, y, frame, f0, df, spec, f_rest,
                vobs=None, rms=None, location=None, mode='w'):
    if vobs is None:
        vobs = 0 * astropy.units.Unit('m/s')
        pass
    
    if rms is None:
        rms = 0 * spec.unit
        pass
    
    array = generate_stage2_array(
        t, 
        x.to('deg').value, 
        y.to('deg').value, 
        f0.to('Hz').value,
        df.to('Hz').value,
        spec.value, 
        vobs.to('m/s').value,
        rms.value,
    )
    
    head = {}
    head['frame'] = frame
    head['spec_unit'] = spec.unit
    head['f_rest'] = f_rest.to('Hz')
    head['location'] = location
    
    bcat.common.bcatio.save(path, array, head, mode)
    return


def load_stage2(path):
    return bcat.common.bcatio.load(path)

