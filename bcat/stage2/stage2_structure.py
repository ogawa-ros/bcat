
import dataclasses
import typing
import numpy
import astropy.units

from bcat.common.common_structure import freq_axis

@dataclasses.dataclass
class stage2_data_structure:
    label: str = ''
    coord: numpy.array = numpy.array([])
    spectrum: numpy.array = numpy.array([])
    rf: freq_axis = freq_axis()
    
    def __post_init__(self):
        if len(self.coord) != self.spectrum.shape[0]:
            errmsg = f'Array length is different between coord and spectrum: ' +\
                     f'({len(self.coord)}, {self.spectrum.shape[0]})'
            raise Exception(errmsg)
        
        pass
        
    def get_rf_axis(self):
        return self.rf.get_axis(self.spectrum.shape[1])

    
