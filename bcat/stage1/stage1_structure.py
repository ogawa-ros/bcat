
import dataclasses
import typing
import numpy
import astropy.units

@dataclasses.dataclass
class stage1_data:
    label: str = ''
    obsmode: typing.List[str] = dataclasses.field(default_factory=list)
    coord: numpy.array = numpy.array([])
    spectrum: numpy.array = numpy.array([])
    rf: freq_axis = freq_axis()
        
    def __post_init__(self):
        if len(self.obsmode) != len(self.coord):
            errmsg = f'Array length is different between obsmode and coord: ' +\
                     f'({len(self.obsmode)}, {len(self.coord)})'
            raise Exception(errmsg)
            
        if len(self.obsmode) != self.spectrum.shape[0]:
            errmsg = f'Array length is different between obsmode and spectrum: ' +\
                     f'({len(self.obsmode)}, {self.spectrum.shape[0]})'
            raise Exception(errmsg)
        
        pass
        
    def get_rf_axis(self):
        return self.rf.get_axis(self.spectrum.shape[1])

    
