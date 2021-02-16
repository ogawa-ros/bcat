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

        hot_mask = numpy.where(numpy.array(self.d.obsmode) == 'HOT')
        off_mask = numpy.where(numpy.array(self.d.obsmode) == 'OFF')
        on_mask = numpy.where(numpy.array(self.d.obsmode) == 'ON')

        spectrum_on = self.d.spectrum[on_mask]
        obstime_on = self.d.coord.obstime[on_mask]

        hot_set = numpy.split(hot_mask[0], numpy.where(numpy.diff(hot_mask[0]) != 1)[0]+1)
        spectrum_hot = numpy.array([numpy.average(self.d.spectrum[hot_set[i]],axis=0) for i in range(len(hot_set)) ])
        obstime_hot  = [self.d.coord.obstime[hot_set[i]][int(len(hot_set[i])/2)] for i in range(len(hot_set))]

        off_set = numpy.split(off_mask[0], numpy.where(numpy.diff(off_mask[0]) != 1)[0]+1)
        spectrum_off = numpy.array([numpy.average(self.d.spectrum[off_set[i]],axis=0) for i in range(len(off_set)) ])
        obstime_off  = [self.d.coord.obstime[off_set[i]][int(len(off_set[i])/2)] for i in range(len(off_set))]

        coord_on = self.d.coord[on_mask]

        obstime_hot  = astropy.time.Time([self.d.coord.obstime[hot_set[i]][int(len(hot_set[i])/2)] for i in range(len(hot_set))])
        obstime_off  = astropy.time.Time([self.d.coord.obstime[off_set[i]][int(len(off_set[i])/2)] for i in range(len(off_set))])

        Tas = []
        for i in range(len(spectrum_on)):
            hot_idx = numpy.abs(obstime_hot-obstime_on[i]).argmin()
            off_idx = numpy.abs(obstime_off-obstime_on[i]).argmin()
            Tas.append(Tamb*(spectrum_on[i]-spectrum_off[off_idx])/(spectrum_hot[hot_idx]-spectrum_off[off_idx]))

        return numpy.array(Tas),coord_on
