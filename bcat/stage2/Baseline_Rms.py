import baseline
import numpy
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

class Baseline_Rms:
    sp_datas = None
    Velocity = None
    mean_spec = None

    fitdata_parts = False

    sigma = None
    maxiters = None
    fitdata_clipped = False

    fit_func = False
    order = None
    chebfit_func = False
    deg = None

    rms = None
    data = None
    baseline_list = None
    rms_list = None
    fit_ch = None

    def __init__(self, sp_datas, Velocity):
        self.Velocity = Velocity.to('km/s').value
        self.sp_datas = numpy.array(sp_datas)
        self.B = baseline.baseline()
        self.B.Velocity = self.Velocity
        self.B.xlim_s = self.Velocity[0]
        self.B.xlim_e = self.Velocity[-1]

        pass

    def mean_sp(self):
        self.mean_spec = numpy.mean(self.sp_datas, axis=0)

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.Velocity, self.mean_spec, label='mean')
        ax.grid(True)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Ta*[K]')
        ax.set_xlim(self.B.xlim_s, self.B.xlim_e)
        ax.legend()
        matplotlib.pyplot.show()

        return

    def input_fitdata_param_parts(self, v_list):
        ##when *v_list
        # x = []
        # for _ in v_list:
        #     x.append(_)

        self.B.data_parts_list = v_list
        self.fitdata_parts = True

        return

    def input_fitdata_param_clipped(self, sigma, maxiters=5):
        self.B.sigma = sigma
        self.B.maxiters = maxiters

        self.fitdata_clipped = True

        return

    def delete_fitdata_param(self):
        self.B.data_parts_list = None
        self.fitdata_parts = False

        self.B.sigma = None
        self.B.maxiters = None
        self.fitdata_clipped = False

        return

    def mk_fit_data(self):
        if self.fitdata_parts:
            if self.fitdata_clipped:
                pass
            else:
                self.B.data_parts()
                self.fit_ch = self.B.fit_parts_ch
                pass
            pass

        elif self.fitdata_clipped:
            if self.fitdata_parts:
                pass
            else:
                self.B.data_clip()
                self.fit_ch = self.B.fit_clip_ch
                pass
            pass

        else:
            pass
        pass

    def select_func(self, fit_func=False, order=None, chebfit_func=False, deg=None):
        self.fit_func = fit_func
        self.order = order
        self.chebfit_func = chebfit_func
        self.deg = deg
        return

    def reset_func(self):
        self.fit_func = False
        self.order = None
        self.chebfit_func = False
        self.deg = None
        return

    def mk_baseline_rms(self):
        if self.fit_func:
            if self.chebfit_func:
                pass
            else:
                self.B.baseline = self.B.baseline_fitting(order=self.order)
                mk_rms_list = []
                for _ in self.fit_ch:
                    mk_rms_list.append(self.B.baseline[_])
                self.rms = numpy.std(mk_rms_list)
                pass
            pass

        elif self.chebfit_func:
            if self.fit_func:
                pass
            else:
                self.B.baseline = self.B.baseline_chebfit(self.deg)
                mk_rms_list = []
                for _ in self.fit_ch:
                    mk_rms_list.append(self.B.baseline[_])
                self.rms = numpy.std(mk_rms_list)

                pass
            pass

        else:
            pass
        pass

    def sample(self):
        self.B.sp_data = self.mean_spec
        self.mk_fit_data()
        print(f'part :{self.fitdata_parts}, clip :{self.fitdata_clipped}')
        print(f'fit v list : {self.B.data_parts_list}')
        print(f'omega :{self.sigma}, maxiters :{self.maxiters}')
        self.B.fit_xy_result()
        self.mk_baseline_rms()
        print(f'x-func :{self.fit_func}, order :{self.order}')
        print(f'cheb func :{self.chebfit_func}, deg :{self.deg}')
        self.B.baseline_graph()
        print(f'rms : {self.rms}')
        return

    def mk_baseline_rms_list(self):
        DATA = []
        data = []
        baseline_list = []
        rms_list = []
        for n in self.sp_datas:
            self.B.sp_data = n
            self.mk_fit_data()
            self.mk_baseline_rms()
            data.append({
            'baseline':self.B.baseline,
            'RMS':self.rms,
            })
            baseline_list.append(self.B.baseline)
            rms_list.append(self.rms)
            self.B.baseline = None
            self.rms = None
            continue

        self.data = data
        self.baseline_list = numpy.array(baseline_list)
        self.rms_list = numpy.array(rms_list)
        DATA.append({
        'baseline_rms':self.data,
        'fit x':self.B.fit_x,
        'fit y':self.B.fit_y,
        'fit ch':self.fit_ch,
        'part':self.fitdata_parts,
        'fit v list':self.B.data_parts_list,
        'clip':self.fitdata_clipped,
        'omega':self.sigma,
        'maxiters':self.maxiters,
        'x-func':self.fit_func,
        'order':self.order,
        'cheb func':self.chebfit_func,
        'deg':self.deg,
        })
        self.DATA = DATA

        return self.DATA, self.baseline_list, self.rms_list, self.fit_ch


    def mk_mean_rms(self):
        rms_list = []
        for i in self.data:
            rms_list.append(i['RMS'])

        return numpy.mean(rms_list)
