import numpy

class stage1_chopper_wheel(object):

    def __init__(self,d1):
        self.d1 = d1
        chopper_wheel()
        pass


    def chopper_wheel(self):

        hot_mask = numpy.where(numpy.array(self.d1.data.obsmode) == 'HOT')
        off_mask = numpy.where(numpy.array(self.d1.data.obsmode) == 'OFF')
        on_mask = numpy.where(numpy.array(self.d1.data.obsmode) == 'ON')

        spectrum_hot = self.d1.data.spectrum[hot_mask]
        spectrum_off = self.d1.data.spectrum[off_mask]
        spectrum_on = self.d1.data.spectrum[on_mask]

        obstime_hot = self.d1.data.coord.obstime[hot_mask]
        obstime_off = self.d1.data.coord.obstime[off_mask]
        obstime_on = self.d1.data.coord.obstime[on_mask]

        coord_on = self.d1.data.coord[on_mask]

        Tas = []
        for i in range(len(spectrum_on)):
            hot_idx = numpy.abs(obstime_hot-obstime_on[i]).argmin()
            off_idx = numpy.abs(obstime_off-obstime_on[i]).argmin()

            Tas.append(Tamb*(spectrum_on[i]-spectrum_off[0])/(spectrum_hot[0]-spectrum_off[0]))

    return numpy.array(Tas),coord_on,obstime_on
