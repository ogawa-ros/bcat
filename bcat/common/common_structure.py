
import dataclasses
import numpy
import astropy.units


dtype_freq = [
    ('cdelt', 'f4'),
    ('crval', 'f4'),
    ('crpix', 'f4'),
]

def freq_axis_array(cdelt, crval, crpix):
    f = numpy.array([cdelt, crval, crpix], dtype='float32').T
    return numpy.frombuffer(f.tobytes(), dtype=dtype_freq)
    

@dataclasses.dataclass
class freq_axis:
    cdelt: astropy.units.Unit = 0 * astropy.units.Hz
    crval: astropy.units.Unit = 0 * astropy.units.Hz
    crpix: int = 0
    
    def __post_init__(self):
        self.cdelt = self.cdelt.to('Hz')
        self.crval = self.crval.to('Hz')
        pass
    
    def get_axis(self, num):
        x = numpy.arange(num)
        return ((x - self.crpix) * self.cdelt + self.crval).to('GHz')


