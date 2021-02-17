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

        mask_hot = numpy.where(numpy.array(self.d.obsmode) == 'HOT')
        mask_off = numpy.where(numpy.array(self.d.obsmode) == 'OFF')
        mask_on = numpy.where(numpy.array(self.d.obsmode) == 'ON')

        set_hot = numpy.split(mask_hot[0], numpy.where(numpy.diff(mask_hot[0]) != 1)[0]+1)
        set_off = numpy.split(mask_off[0], numpy.where(numpy.diff(mask_off[0]) != 1)[0]+1)

        spectrum_hot = numpy.array([numpy.average(self.d.spectrum[set_hot[i]],axis=0) for i in range(len(set_hot)) ])
        spectrum_off = numpy.array([numpy.average(self.d.spectrum[set_off[i]],axis=0) for i in range(len(set_off)) ])
        spectrum_on = self.d.spectrum[mask_on]

        obstime_hot  = astropy.time.Time([self.d.coord.obstime[set_hot[i]][int(len(set_hot[i])/2)] for i in range(len(set_hot))])
        obstime_off  = astropy.time.Time([self.d.coord.obstime[off_set[i]][int(len(set_off[i])/2)] for i in range(len(set_off))])
        obstime_on = self.d.coord.obstime[mask_on]

        coord_on = self.d.coord[mask_on]


        _Tas = []
        for i in range(len(spectrum_on)):
            idx_hot = numpy.abs(obstime_hot-obstime_on[i]).argmin()
            idx_off = numpy.abs(obstime_off-obstime_on[i]).argmin()
            Tas.append(Tamb*(spectrum_on[i]-spectrum_off[idx_off])/(spectrum_hot[idx_hot]-spectrum_off[idx_off]))
        Tas = numpy.array(_Tas)

        return Tas,coord_on
