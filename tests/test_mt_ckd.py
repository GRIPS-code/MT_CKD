from os.path import dirname, join, realpath

import matplotlib.pyplot as plt
from netCDF4 import Dataset
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
    scriptdir = dirname(realpath(__file__))
    path = join(scriptdir, "../mt_ckd/mt-ckd.nc")
    temperature = 296. # [K].
    pressure = 1013. # [mb].
    volume_mixing_ratio = {"Ar": 0.009,
                           "CO2": 345.e-6,
                           "H2O": 0.01,
                           "N2": 0.78,
                           "O2": 0.21,
                           "O3": 3.75e-6,} #[mol mol-1].
    grid = arange(0., 10002., 2)
    extinction = zeros(grid.size)
#   path_length = 1. # [cm].
    continua = {"CO2": CarbonDioxideContinuum,
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
    test_results = {"H2O_s": "h2o-self", "H2O_f": "h2o-foreign", "CO2": "co2",
                    "O3": "o3", "O2": "o2", "N2": "n2"}

    for key, value in test_results.items():
        with Dataset(join(scriptdir, "{}.nc".format(value)), "r") as dataset:
            test_grid = dataset.variables["wavenumber"][:]
            test_data = dataset.variables["absrb"][:]
            plt.plot(test_grid, test_data, color="black", label="original")

        if key in ["CO2", "H2O_f", "H2O_s"]:
            keys = [key,]
        else:
            keys = [x for x in continua.keys() if x.startswith(key)]
        for i in keys:
            x = continua[i](path)
            xp = x.grid()
            fp = x.spectra(temperature, pressure, volume_mixing_ratio)
#           fp2 = interp(grid, xp, fp, left=0., right=0.)
            plt.plot(xp, fp, color="red", linestyle="--", label="python")

        plt.yscale("log")
        plt.legend()
        plt.title(key)
        plt.show()
