from os.path import dirname, join, realpath

import matplotlib.pyplot as plt

from mt_ckd import CarbonDioxideContinuum, \
                   NitrogenCIAPureRotationContinuum, NitrogenCIAFundamentalContinuum, \
                   NitrogenCIAFirstOvertoneContinuum, \
                   OxygenCIAFundamentalContinuum, OxygenCIANIRContinuum, \
                   OxygenCIANIR2Continuum, OxygenCIANIR3Continuum, OxygenVisibleContinuum, \
                   OxygenHerzbergContinuum, OxygenUVContinuum, \
                   OzoneChappuisWulfContinuum, OzoneHartleyHugginsContinuum, OzoneUVContinuum, \
                   WaterVaporForeignContinuum, WaterVaporSelfContinuum


if __name__ == "__main__":
    #User inputs.
    path = join(dirname(realpath(__file__)), "../mt_ckd/mt-ckd.nc")
    temperature = 260. # [K].
    pressure = 1013. # [mb].
    volume_mixing_ratio = {"Ar": 0.009,
                           "CO2": 345.e-6,
                           "H2O": 0.01,
                           "N2": 0.78,
                           "O2": 0.21,
                           "O3": 3.75e-6,} #[mol mol-1].
#   path_length = 1. # [cm].

    #Carbon dioxide.
    co2 = CarbonDioxideContinuum(path)
    plt.plot(co2.grid(), co2.spectra(temperature, pressure, volume_mixing_ratio),
             color="red", label="CO2")

    #Water vapor.
    h2o_self = WaterVaporSelfContinuum(path)
    plt.plot(h2o_self.grid(), h2o_self.spectra(temperature, pressure, volume_mixing_ratio),
             color="orange", label="H2O self")
    h2o_foreign = WaterVaporForeignContinuum(path)
    plt.plot(h2o_foreign.grid(), h2o_foreign.spectra(temperature, pressure, volume_mixing_ratio),
             color="yellow", label="H2O foreign")

    #Ozone.
    o3_cw = OzoneChappuisWulfContinuum(path)
    plt.plot(o3_cw.grid(), o3_cw.spectra(temperature, pressure, volume_mixing_ratio),
             color="blue", label="O3_cw")
    o3_hh = OzoneHartleyHugginsContinuum(path)
    plt.plot(o3_hh.grid(), o3_hh.spectra(temperature, pressure, volume_mixing_ratio),
             color="blue", label="O3_hh")
    o3_uv = OzoneUVContinuum(path)
    plt.plot(o3_uv.grid(), o3_uv.spectra(temperature, pressure, volume_mixing_ratio),
             color="blue", label="O3_uv")

    #Oxygen
    o2_f = OxygenCIAFundamentalContinuum(path)
    plt.plot(o2_f.grid(), o2_f.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_f")
    o2_cia_nir = OxygenCIANIRContinuum(path)
    plt.plot(o2_cia_nir.grid(), o2_cia_nir.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_cia_nir")
    o2_cia_nir2 = OxygenCIANIR2Continuum()
    plt.plot(o2_cia_nir2.grid(), o2_cia_nir2.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_cia_nir2")
    o2_cia_nir3 = OxygenCIANIR3Continuum(path)
    plt.plot(o2_cia_nir3.grid(), o2_cia_nir3.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_cia_nir3")
    o2_cia_vis = OxygenVisibleContinuum(path)
    plt.plot(o2_cia_vis.grid(), o2_cia_vis.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_cia_vis")
    o2_herz = OxygenHerzbergContinuum()
    plt.plot(o2_herz.grid(), o2_herz.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_herz")
    o2_uv = OxygenUVContinuum(path)
    plt.plot(o2_uv.grid(), o2_uv.spectra(temperature, pressure, volume_mixing_ratio),
             color="green", label="O2_uv")

    #Nitrogen
    n2_r = NitrogenCIAPureRotationContinuum(path)
    plt.plot(n2_r.grid(), n2_r.spectra(temperature, pressure, volume_mixing_ratio),
             color="black", label="N2_r")
    n2_f = NitrogenCIAFundamentalContinuum(path)
    plt.plot(n2_f.grid(), n2_f.spectra(temperature, pressure, volume_mixing_ratio),
             color="black", label="N2_f")
    n2_o = NitrogenCIAFirstOvertoneContinuum(path)
    plt.plot(n2_o.grid(), n2_o.spectra(temperature, pressure, volume_mixing_ratio),
             color="black", label="N2_o")

    plt.yscale("log")
    plt.legend()
    plt.show()
