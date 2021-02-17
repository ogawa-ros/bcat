import numpy
import astropy
import time

class stage1_chopper_wheel(object):

    def __init__(self,d):
        self.d = d
        pass


    def chopper_wheel(self):
        _Tamb = self.d.coord.temperature
        Tamb = _Tamb.to(astropy.units.K, equivalencies=astropy.units.temperature())

        mask_hot = numpy.where(numpy.array(self.d.obsmode) == 'HOT')[0]
        mask_off = numpy.where(numpy.array(self.d.obsmode) == 'OFF')[0]
        mask_on = numpy.where(numpy.array(self.d.obsmode) == 'ON')[0]

        set_hot = numpy.split(mask_hot, numpy.where(numpy.diff(mask_hot) != 1)[0]+1)
        set_off = numpy.split(mask_off, numpy.where(numpy.diff(mask_off) != 1)[0]+1)
        #それぞれのsetの最後の値のindex
        _idx_hot = mask_hot[numpy.where(numpy.diff(mask_hot) != 1)]
        _idx_off = mask_off[numpy.where(numpy.diff(mask_off) != 1)]
        _idx_hot = numpy.append(_idx_hot,mask_hot[-1])
        _idx_off = numpy.append(_idx_off,mask_off[-1])

        spectrum_hot = numpy.array([numpy.average(self.d.spectrum[set_hot[i]],axis=0) for i in range(len(set_hot)) ])
        spectrum_off = numpy.array([numpy.average(self.d.spectrum[set_off[i]],axis=0) for i in range(len(set_off)) ])
        spectrum_on = self.d.spectrum[mask_on]

        coord_on = self.d.coord[mask_on]


        _Tas = []
        for i in range(len(mask_on)):
            idx_hot = numpy.abs(_idx_hot-mask_on[i]).argmin()
            idx_off = numpy.abs(_idx_off-mask_on[i]).argmin()
            _Tas.append(Tamb*(spectrum_on[i]-spectrum_off[idx_off])/(spectrum_hot[idx_hot]-spectrum_off[idx_off]))
        Tas = numpy.array(_Tas)

        return Tas,coord_on
