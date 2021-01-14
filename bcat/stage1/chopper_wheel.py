import numpy
import astropy

class stage1_chopper_wheel(object):

    def __init__(self,d1):
        self.d1 = d1
        pass


    def chopper_wheel(self):
        _Tamb = self.d1.coord.temperature
        Tamb = _Tamb.to(astropy.units.K, equivalencies=astropy.units.temperature())

        hot_mask = numpy.where(numpy.array(self.d1.obsmode) == 'HOT')
        off_mask = numpy.where(numpy.array(self.d1.obsmode) == 'OFF')
        on_mask = numpy.where(numpy.array(self.d1.obsmode) == 'ON')

        spectrum_hot = self.d1.spectrum[hot_mask]
        spectrum_off = self.d1.spectrum[off_mask]
        spectrum_on = self.d1.spectrum[on_mask]

        obstime_hot = self.d1.coord.obstime[hot_mask]
        obstime_off = self.d1.coord.obstime[off_mask]
        obstime_on = self.d1.coord.obstime[on_mask]

        coord_on = self.d1.coord[on_mask]

        Tas = []
        for i in range(len(spectrum_on)):
            hot_idx = numpy.abs(obstime_hot-obstime_on[i]).argmin()
            off_idx = numpy.abs(obstime_off-obstime_on[i]).argmin()

            Tas.append(Tamb*(spectrum_on[i]-spectrum_off[off_idx])/(spectrum_hot[hot_idx]-spectrum_off[off_idx]))

        return numpy.array(Tas),coord_on
