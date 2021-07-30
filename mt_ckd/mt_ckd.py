from netCDF4 import Dataset
from numpy import copy, power, zeros


class Spectrum(object):
    def __init__(self, path, name):
        with Dataset(path, "r") as dataset:
            v = dataset.variables[name]
            self.data = copy(v[:])
            self.grid = {x: v.getncattr("wavenumber_{}".format(x)) for x in
                         ["lower_bound", "upper_bound", "resolution"]}
#           self.units = v.getncattr("units")

    def wavenumbers(self):
        return [self.grid["lower_bound"] + i*self.grid["resolution"]
                for i in range(self.data.size)]


class Continuum(object):
    def __init__(self, path):
        raise NotImplementedError("You must override this class.")

    def spectra(self):
        raise NotImplementedError("You must override this method.")


class CarbonDioxideContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "bfco2")
        self.t_correction = [1.44e-01, 3.61e-01, 5.71e-01, 7.63e-01, 8.95e-01,
                             9.33e-01, 8.75e-01, 7.30e-01, 5.47e-01, 3.79e-01,
                             2.55e-01, 1.78e-01, 1.34e-01, 1.07e-01, 9.06e-02,
                             7.83e-02, 6.83e-02, 6.00e-02, 5.30e-02, 4.72e-02,
                             4.24e-02, 3.83e-02, 3.50e-02, 3.23e-02, 3.01e-02]
        self.xfac_co2 = Spectrum(path, "x_factor_co2")

    def spectra(self, temperature):
        scale = zeros(self.data.data.size)
        t_factor = zeros(self.data.data.size)
        for i in range(self.data.data.size):
            vj = self.data.grid["lower_bound"] + i*self.data.grid["resolution"]
            if vj >= 2000. and vj <= 2998.:
                jfac = int((vj - 1998.)/2. + 0.00001) - 1
                scale[i] = self.xfac_co2.data[jfac]
            else:
                scale[i] = 1.
            if vj >= 2386. and vj <= 2424.:
                jfac = int((vj - 2386.)/2.)
                t_factor[i] = power(temperature/246., self.t_correction[jfac])
            else:
                t_factor[i] = 1.
        return scale[:]*t_factor[:]*self.data.data[:]

    def grid(self):
        return self.data.wavenumbers()


class OzoneContinuum(Continuum):
    def __init__(self, path):
        self.data = None

    def spectra(self, temperature):
        pass



class WaterVaporForeignContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "bfh2o")
        self.xfac_rhu = [0.7620, 0.7840, 0.7820, 0.7840, 0.7620, 0.7410, 0.7970,
                         0.7810, 0.8330, 0.8500, 0.8330, 0.7810, 0.7540, 0.8180,
                         0.9140, 0.9980, 0.9830, 0.9330, 0.8850,
                         0.8420, 0.8070, 0.8000, 0.8010, 0.8100,
                         0.8090, 0.8320, 0.8180, 0.7970, 0.8240,
                         0.8640, 0.8830, 0.8830, 0.8470, 0.8380,
                         0.8660, 0.9410, 1.0400, 1.0680, 1.1410,
                         1.0800, 1.0340, 1.1550, 1.0990, 1.0270,
                         0.9500, 0.8950, 0.8150, 0.7830, 0.7700,
                         0.7000, 0.7650, 0.7750, 0.8500, 0.9000,
                         0.9050, 0.9540, 1.0200, 1.0200, 1.0250,
                         1.0200, 1.1000, 1.1250, 1.1200, 1.1110,
                         1.1370, 1.1600, 1.1490, 1.1070, 1.0640,
                         1.0450]

    def spectra(self):
        scale = zeros(self.data.data.size)
        for i in range(self.data.data.size):
            vj = self.data.grid["lower_bound"] + i*self.data.grid["resolution"]
            if vj <= 600.:
                jfac = int((vj + 10.)/10. +  0.00001) + 1
                scale[i] = self.xfac_rhu[jfac]
            else:
                vdelsq1 = (vj - 255.67)*(vj - 255.67)
                vdelmsq1 = (vj + 255.67)*(vj + 255.67)
                vf1 = power((vj - 255.67)/57.83, 8)
                vmf1 = power((vj + 255.67)/57.83, 8)
                vf2 = power(vj/630., 8)
                scale[i] = 1. + (0.06 + -0.42*(((240.*240.)/(vdelsq1 + 240.*240. + vf1)) +
                           ((240.*240.)/(vdelmsq1 + 240.*240. + vmf1))))/(1. + 0.3*vf2)
        return 1.e-20*scale[:]*self.data.data[:]

    def grid(self):
        return self.data.wavenumbers()


class WaterVaporSelfContinuum(Continuum):
    def __init__(self, path):
        self.data = {296: Spectrum(path, "bs296"),
                     260: Spectrum(path, "bs260")}

    def spectra(self, temperature):
        t_factor = (temperature - 296.)/(260. - 296.)
        return 1.e-20*self.data[296].data[:] * \
               power(self.data[260].data[:]/self.data[296].data[:], t_factor)

    def grid(self):
        return self.data[296].wavenumbers()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    h2o_self = WaterVaporSelfContinuum("mt-ckd.nc")
    plt.plot(h2o_self.grid(), h2o_self.spectra(300.), label="H2O self")
    h2o_foreign = WaterVaporForeignContinuum("mt-ckd.nc")
    plt.plot(h2o_foreign.grid(), h2o_foreign.spectra(), label="H2O foreign")
    co2 = CarbonDioxideContinuum("mt-ckd.nc")
    plt.plot(co2.grid(), co2.spectra(300.), label="CO2")
    plt.yscale("log")
    plt.legend()
    plt.show()
