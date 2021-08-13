from os.path import dirname, join, realpath

import matplotlib.pyplot as plt
from numpy import arange, interp, zeros

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
    temperature = 296. # [K].
    pressure = 500. # [mb].
    volume_mixing_ratio = {"Ar": 0.009,
                           "CO2": 345.e-6,
                           "H2O": 0.01,
                           "N2": 0.78,
                           "O2": 0.21,
                           "O3": 3.75e-6,} #[mol mol-1].
    grid = arange(0., 10002., 2)
    extinction = zeros(grid.size)
#   path_length = 1. # [cm].
    continuua = {"CO2": CarbonDioxideContinuum,
                 "N2_rot": NitrogenCIAPureRotationContinuum,
                 "N2_fund": NitrogenCIAFundamentalContinuum,
                 "N2_over": NitrogenCIAFirstOvertoneContinuum,
                 "O2_fund": OxygenCIAFundamentalContinuum,
                 "O2_nir1": OxygenCIANIRContinuum,
                 "O2_nir2": OxygenCIANIR2Continuum,
                 "O2_nir3": OxygenCIANIR3Continuum,
                 "O2_vis": OxygenVisibleContinuum,
                 "O2_herz": OxygenHerzbergContinuum,
                 "O2_uv": OxygenUVContinuum,
                 "O3_cw": OzoneChappuisWulfContinuum,
                 "O3_hh": OzoneHartleyHugginsContinuum,
                 "O3_uv": OzoneUVContinuum,
                 "H2O_f": WaterVaporForeignContinuum,
                 "H2O_s": WaterVaporSelfContinuum}

    for label, continuum in continuua.items():
        #Plot the component.
        x = continuum(path)
        xp = x.grid()
        fp = x.spectra(temperature, pressure, volume_mixing_ratio)
        plt.plot(xp, fp,label=label)

        #Interpolate to input grid.
        fp2 = interp(grid, xp, fp, left=0., right=0.)
        extinction = extinction[:] + fp2[:]
        plt.plot(grid, extinction)

    plt.plot(grid, extinction, label="total")
    plt.yscale("log")
    plt.legend()
    plt.show()
