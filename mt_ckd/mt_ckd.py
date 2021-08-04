from netCDF4 import Dataset
from numpy import arange, asarray, copy, exp, log, power, zeros


class Spectrum(object):
    """Helper class that reads data from a variable in the input dataset.

    Attributes:
        path: Path to the netcdf dataset.
        grid: Dictionary describing the wavenumber grid.
    """
    def __init__(self, path, name):
        """Reads the data from a variable in the input dataset.

        Args:
            path: Path to the netcdf dataset.
            name: Name of the variable in the dataset.
        """
        with Dataset(path, "r") as dataset:
            v = dataset.variables[name]
            self.data = copy(v[:])
            self.grid = {x: v.getncattr("wavenumber_{}".format(x)) for x in
                         ["lower_bound", "upper_bound", "resolution"]}
#           self.units = v.getncattr("units")

    def wavenumbers(self):
        """Cacluates the wavenumber grid [cm-1] for the variable.

        Returns:
            A 1d numpy array containing the wavenumber grid [cm-1].
        """
        return asarray([self.grid["lower_bound"] + i*self.grid["resolution"]
                        for i in range(self.data.size)])


class Continuum(object):
    """Abstract class for gridded continuum coefficients."""
    def __init__(self, path):
        """Reads in the necessary data from an input dataset.

        Args:
            path: Path to the netcdf dataset.
        """
        raise NotImplementedError("You must override this class.")

    def spectra(self, temperature):
        """Calculates the continuum coefficients at an input temperature.

        Args:
            temperature: Temperaure [K].

        Returns:
            An array of continuum coefficients.
        """
        raise NotImplementedError("You must override this method.")

    def grid(self):
        """Calculates the wavenumber grid [cm-1].

        Returns:
            A 1d numpy array containing the wavenumber grid [cm-1].
        """
        raise NotImplementedError("You must override this method.")


class CarbonDioxideContinuum(Continuum):
    """Carbon dioxide continuum coefficients.

    Attributes:
        data: Spectrum object containing data read from an input dataset.
        t_correction: Temperature correction coefficients.
        xfac_co2: Spectrum object of chi-factors?
    """
    def __init__(self, path):
        self.data = Spectrum(path, "bfco2")
        self.t_correction = [1.44e-01, 3.61e-01, 5.71e-01, 7.63e-01, 8.95e-01,
                             9.33e-01, 8.75e-01, 7.30e-01, 5.47e-01, 3.79e-01,
                             2.55e-01, 1.78e-01, 1.34e-01, 1.07e-01, 9.06e-02,
                             7.83e-02, 6.83e-02, 6.00e-02, 5.30e-02, 4.72e-02,
                             4.24e-02, 3.83e-02, 3.50e-02, 3.23e-02, 3.01e-02]
        self.xfac_co2 = Spectrum(path, "x_factor_co2")

    def spectra(self, temperature, pressure):
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
        return 1.e-20*(pressure/1013.)*(296./temperature)*scale[:]* \
               t_factor[:]*self.data.data[:]

    def grid(self):
        return self.data.wavenumbers()


class NitrogenCIAPureRotationContinuum(Continuum):
    def __init__(self, path):
        self.data = {296: [Spectrum(path, "ct_296"), Spectrum(path, "sf_296")],
                     220: [Spectrum(path, "ct_220"), Spectrum(path, "sf_220")]}

    def spectra(self, temperature, pressure, vmr):
        factor = (temperature - 296.)/(220. - 296.)
        c = self.data[296][0].data[:]*power(self.data[220][0].data[:]/self.data[296][0].data[:],
                                            factor)
        s = self.data[296][1].data[:]*power(self.data[220][1].data[:]/self.data[296][1].data[:],
                                            factor)
        fo2 = (s[:] - 1.)*0.79/0.21
        return (1./2.68675e19)*(pressure/1013.)*(273./temperature)*c[:]* \
               (vmr["N2"] + fo2[:]*vmr["O2"] + vmr["H2O"])

    def grid(self):
        return self.data[296][0].wavenumbers()


class NitrogenCIAFundamentalContinuum(Continuum):
    def __init__(self, path):
        self.data = [Spectrum(path, "xn2_272"), Spectrum(path, "xn2_228"),
                     Spectrum(path, "a_h2o")]

    def spectra(self, temperature, pressure, vmr):
        xtfac = ((1./temperature) - (1./272.))/((1./228.) - (1./272.))
        xt_lin = (temperature - 272.)/(228. - 272.)
        factor = (1./2.68675e19)*(pressure/1013.)*(273./temperature)
        c0 = zeros(self.data[0].data.size)
        for i in range(self.data[0].data.size):
            if self.data[0].data[i] > 0. and self.data[1].data[i] > 0.:
                c0[i] = self.data[0].data[i]*power(self.data[1].data[i]/self.data[0].data[i],
                                                   xtfac)
            else:
                c0[i] = self.data[0].data[i] + (self.data[1].data[i] - self.data[0].data[i])* \
                        xt_lin
        c1 = (1.294 - 0.4545*temperature/296.)*c0[:]
        c2 = (9./7.)*self.data[2].data[:]*c0[:]
        return factor*(vmr["N2"]*c0[:] + vmr["O2"]*c1[:] + vmr["H2O"]*c2[:])/self.grid()[:]

    def grid(self):
        return self.data[0].wavenumbers()


class NitrogenCIAFirstOvertoneContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "xn2")

    def spectra(self, temperature, pressure, vmr):
        factor = (1./2.68675e19)*(pressure/1013.)*(273./temperature)* \
                 (vmr["N2"] + vmr["O2"] + vmr["H2O"])
        return factor*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OzoneChappuisWulfContinuum(Continuum):
    """Ozone continuum in the Chappuis and Wulf band.

    Attributes:
        data: List of Spectrum objects containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = [Spectrum(path, "x_o2"), Spectrum(path, "y_o2"), Spectrum(path, "z_o2")]

    def spectra(self, temperature):
        dt = temperature - 273.15
        return 1.e-20*(self.data[0].data[:] + self.data[1].data[:]*dt +
                       self.data[2].data[:]*dt*dt)/self.grid()[:]

    def grid(self):
        return self.data[0].wavenumbers()


class OzoneHartleyHugginsContinuum(Continuum):
    """Ozone Hartly-Huggins continuum cros sections.

    Attributes:
        data: List of Spectrum objects containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = [Spectrum(path, "o3_hh0"), Spectrum(path, "o3_hh1"),
                     Spectrum(path, "o3_hh2")]

    def spectra(self, temperature):
        dt = temperature - 273.15
        return 1.e-20*(self.data[0].data[:]/self.grid()[:])* \
               (1. + self.data[1].data[:]*dt + self.data[2].data[:]*dt*dt)

    def grid(self):
        return self.data[0].wavenumbers()


class OzoneUVContinuum(Continuum):
    """Ozone ultra-violet continuum coefficients.

    Attributes:
        data: A Spectrum object containing data read from an input dataset.
    """
    def __init__(self, path):
        self.data = Spectrum(path, "o3_huv")

    def spectra(self):
        return self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenCIAFundamentalContinuum(Continuum):
    def __init__(self, path):
        self.data = [Spectrum(path, "o2_f"), Spectrum(path, "o2_t")]

    def spectra(self, temperature, pressure):
        """Pressure in mb?"""
        xktfac = (1./296.) - (1./temperature)
        return 1.e-20*(pressure/1013.)*(273./temperature)*(1.e20/2.68675e19)* \
               self.data[0].data[:]*exp(self.data[1].data[:]*xktfac)/ \
               self.data[0].wavenumbers()[:]

    def grid(self):
        return self.data[0].wavenumbers()


class OxygenCIANIRContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_inf1")

    def spectra(self, temperature, pressure, vmr):
        factor = (1./2.68675e19)*(pressure/1013.)*(273./temperature)* \
                 ((1./0.446)*vmr["O2"] + (0.3/0.446)*vmr["N2"] + vmr["H2O"])
        return factor*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenCIANIR2Continuum(Continuum):
    def __init__(self, path=None):
        self._grid = arange(9100., 11002., 2.)
        self.data = zeros(self._grid.size)
        hw1 = 58.96
        hw2 = 45.04
        for i in range(self._grid.size):
            dv1 = self._grid[i] - 9375.
            dv2 = self._grid[i] - 9439.
            damp1 = exp(dv1/176.1) if dv1 < 0. else 1.
            damp2 = exp(dv2/176.1) if dv2 < 0. else 1.
            o2inf = 0.31831*(((1.166e-04*damp1/hw1)/(1. + (dv1/hw1)*(dv1/hw1))) +
                             ((3.086e-05*damp2/hw2)/(1. + (dv2/hw2)*(dv2/hw2))))*1.054
            self.data[i] = o2inf/self._grid[i]

    def spectra(self, temperature, pressure):
        adjwo2 = 1.e-20*(pressure/1013.)*(296./temperature)
        return self.data[:]*adjwo2

    def grid(self):
        return self._grid[:]


class OxygenCIANIR3Continuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_inf3")

    def spectra(self, temperature, pressure):
        factor = (1./2.68675e+19)*(pressure/1013.)*(273./temperature)
        return factor*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenVisibleContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_invis")

    def spectra(self, temperature, pressure):
        factor = 1./((2.68675e19*1.e-20*(55.*273./296.)*(55.*273./296.))*89.5)
        adjwo2 = 1.e-20*(pressure/1013.)*(273./temperature)
        return adjwo2*factor*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class OxygenHerzbergContinuum(Continuum):
    def __init__(self, path=None):
        self._grid = arange(36000., 100010., 10.)
        self.data = zeros(self._grid.size)
        for i in range(self._grid.size):
            if self._grid[i] <= 36000.:
                self.data[i] = 0.
            else:
                corr = ((40000. - self._grid[i])/4000.)*7.917e-7 \
                       if self._grid[i] <= 40000. else 0.
                yratio = self._grid[i]/48811.0
                self.data[i] = 6.884e-4*yratio*exp(-69.738*power(log(yratio), 2)) - corr

    def spectra(self, temperature, pressure):
        factor = 1. + 0.83*(pressure/1013.)*(273.16/temperature)
        return 1.e-20*factor*self.data[:]/self.grid()[:]

    def grid(self):
        return self._grid[:]


class OxygenUVContinuum(Continuum):
    def __init__(self, path):
        self.data = Spectrum(path, "o2_infuv")

    def spectra(self):
        return 1.e-20*self.data.data[:]/self.grid()[:]

    def grid(self):
        return self.data.wavenumbers()


class WaterVaporForeignContinuum(Continuum):
    """Water vapor foreign continuum coefficients.

    Attributes:
        data: Spectrum object containing data read from an input dataset.
        xfac_rhu: List of ???
    """
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
    """Water vapor self continuum coefficients.

    Attributes:
        data: Dictionary that maps temperatures (keys) to Spectrum objects containing
              data read from an input dataset (values).
    """
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

    #User inputs.
    path = "mt-ckd.nc"
    temperature = 260. # [K].
    pressure = 1013. # [mb].
    volume_mixing_ratio = {"Ar": 0.009,
                           "CO2": 345.e-6,
                           "H2O": 0.01,
                           "N2": 0.78,
                           "O2": 0.21,
                           "O3": 3.75e-6,}
    path_length = 1. # [cm].

    #MT-CKD required variables.
    w = {"dry": (pressure/1013.)*(273./temperature)*path_length*
                (1. - volume_mixing_ratio["H2O"]),}
    for key, value in volume_mixing_ratio.items():
        w[key] = w["dry"]*value
    w["broad"] = w["N2"] + w["Ar"]
    w["total"] = sum([w[x] for x in volume_mixing_ratio.keys()])
    vmr = {"H2O": w["H2O"]/w["total"], "O2": w["O2"]/w["total"]}
    vmr["N2"] = 1. - vmr["H2O"] - vmr["O2"]

    #Carbon dioxide.
    co2 = CarbonDioxideContinuum(path)
    plt.plot(co2.grid(), w["CO2"]*co2.spectra(temperature, pressure),
             color="red", label="CO2")

    #Water vapor.
    h2o_self = WaterVaporSelfContinuum(path)
    plt.plot(h2o_self.grid(), h2o_self.spectra(temperature), color="orange", label="H2O self")
    h2o_foreign = WaterVaporForeignContinuum(path)
    plt.plot(h2o_foreign.grid(), h2o_foreign.spectra(), color="yellow", label="H2O foreign")

    #Ozone.
    o3_cw = OzoneChappuisWulfContinuum(path)
    plt.plot(o3_cw.grid(), w["O3"]*o3_cw.spectra(temperature), color="blue",
             label="O3_cw")
    o3_hh = OzoneHartleyHugginsContinuum(path)
    plt.plot(o3_hh.grid(), w["O3"]*o3_hh.spectra(temperature), color="blue", label="O3_hh")
    o3_uv = OzoneUVContinuum(path)
    plt.plot(o3_uv.grid(), w["O3"]*o3_uv.spectra(), color="blue", label="O3_uv")

    #Oxygen
    o2_f = OxygenCIAFundamentalContinuum(path)
    plt.plot(o2_f.grid(), w["O2"]*o2_f.spectra(temperature, pressure), color="green",
             label="O2_f")
    o2_cia_nir = OxygenCIANIRContinuum(path)
    plt.plot(o2_cia_nir.grid(), w["O2"]*o2_cia_nir.spectra(temperature, pressure, vmr),
             color="green", label="O2_cia_nir")
    o2_cia_nir2 = OxygenCIANIR2Continuum()
    plt.plot(o2_cia_nir2.grid(),
             w["O2"]*(w["O2"]/w["total"])*(1./0.209)*o2_cia_nir2.spectra(temperature, pressure),
             color="green", label="O2_cia_nir2")
    o2_cia_nir3 = OxygenCIANIR3Continuum(path)
    plt.plot(o2_cia_nir3.grid(), w["O2"]*o2_cia_nir3.spectra(temperature, pressure),
             color="green", label="O2_cia_nir3")
    o2_cia_vis = OxygenVisibleContinuum(path)
    plt.plot(o2_cia_vis.grid(),
             w["O2"]*(w["O2"]/w["total"])*o2_cia_vis.spectra(temperature, pressure),
             color="green", label="O2_cia_vis")
    o2_herz = OxygenHerzbergContinuum()
    plt.plot(o2_herz.grid(), w["O2"]*o2_herz.spectra(temperature, pressure),
             color="green", label="O2_herz")
    o2_uv = OxygenUVContinuum(path)
    plt.plot(o2_uv.grid(), w["O2"]*o2_uv.spectra(), color="green", label="O2_uv")

    #Nitrogen
    n2_r = NitrogenCIAPureRotationContinuum(path)
    plt.plot(n2_r.grid(),
             vmr["N2"]*w["total"]*n2_r.spectra(temperature, pressure, vmr),
             color="black", label="N2_r")
    n2_f = NitrogenCIAFundamentalContinuum(path)
    plt.plot(n2_f.grid(),
             vmr["N2"]*w["total"]*n2_f.spectra(temperature, pressure, vmr),
             color="black", label="N2_f")
    n2_o = NitrogenCIAFirstOvertoneContinuum(path)
    plt.plot(n2_o.grid(),
             vmr["N2"]*w["total"]*n2_o.spectra(temperature, pressure, vmr),
             color="black", label="N2_o")

    plt.yscale("log")
    plt.legend()
    plt.show()
