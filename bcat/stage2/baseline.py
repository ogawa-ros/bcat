import numpy
from astropy.stats import sigma_clip, mad_std
import matplotlib.pyplot
#%matplotlib inline
matplotlib.rcParams['image.origin'] = 'lower'
matplotlib.rcParams['image.interpolation'] = 'none'
matplotlib.rcParams['image.cmap'] = 'inferno'
matplotlib.rcParams['font.family'] = 'Arial,freesans'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.left'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['ytick.direction'] = 'in'

class baseline:
    Velocity = None
    sp_data = None
    fit_x = None
    fit_y = None
    baseline = None

    data_parts_list = None
    fit_parts_ch = None

    sigma = None
    maxiters = 5
    fit_clip_ch = None

    xlim_s = None
    xlim_e = None

    #def __init__(self, sp_data, Velocity):
    #    self.Velocity = Velocity
        #self.sp_data = numpy.array(sp_data)

    #    pass

    def data_parts(self):
        V = self.Velocity
        x = []
        y = []
        ch = []
        for _ in self.data_parts_list:
            fitnum1 = numpy.argmin(abs(V-_[0]))
            fitnum2 = numpy.argmin(abs(V-_[1]))
            x.append(V[fitnum1:fitnum2+1])
            y.append(self.sp_data[fitnum1:fitnum2+1])
            ch += list(range(len(V))[fitnum1:fitnum2+1])

        self.fit_x = numpy.concatenate(x)
        self.fit_y = numpy.concatenate(y)
        self.fit_parts_ch = numpy.array(ch)

        return self.fit_x, self.fit_y

    def data_clip(self):
        f = self.Velocity
        filtered_data = sigma_clip(self.sp_data,
                                   sigma=self.sigma,
                                   maxiters=self.maxiters,
                                   stdfunc=mad_std,
                                   masked=True)
        self.fit_x = f[~(filtered_data.mask)]
        self.fit_y = self.sp_data[~(filtered_data.mask)]
        self.fit_clip_ch = numpy.array(range(len(f)))[~(filtered_data.mask)]

        return self.fit_x, self.fit_y


    def fitting_func3(self, x, a, b, c, d):
        return a*x**3 + b*x**2 + c*x**1 + d

    def fitting_func2(self, x, a, b, c):
        return a*x**2 + b*x**1 + c

    def fitting_func1(self, x, a, b):
        return a*x + b

    def baseline_fitting(self, order=None):
        param = numpy.polyfit(self.fit_x, self.fit_y, order)
        if order == 1:
            self.baseline = self.sp_data - self.fitting_func1(self.Velocity, *param)
        elif order == 2:
            self.baseline = self.sp_data - self.fitting_func2(self.Velocity, *param)
        elif order == 3:
            self.baseline = self.sp_data - self.fitting_func3(self.Velocity, *param)
        else:
            pass

        return self.baseline


    def fitting_func_cheb(self, array_x, array_y, X, deg):
        c = numpy.polynomial.chebyshev.chebfit(array_x, array_y, deg=deg)

        return numpy.polynomial.chebyshev.chebval(X, c)

    def baseline_chebfit(self, deg):
        cheb = self.fitting_func_cheb(self.fit_x, self.fit_y, self.Velocity, deg)
        self.baseline = self.sp_data - cheb

        return self.baseline

    def fit_xy_result(self):
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.Velocity, self.sp_data, label='raw data')
        ax.plot(self.fit_x, self.fit_y, '.', label='fit data')
        ax.grid(True)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Ta*[K]')
        ax.legend()
        matplotlib.pyplot.show()

        return

    def baseline_graph(self):
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.Velocity, self.sp_data, label='raw data')
        ax.plot(self.Velocity, self.baseline, label='baseline')
        ax.grid(True)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Ta*[K]')
        ax.set_xlim(self.xlim_s, self.xlim_e)
        ax.legend()
        matplotlib.pyplot.show()

        return
