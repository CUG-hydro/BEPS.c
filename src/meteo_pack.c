/// @file meteo_pack.c
/// @brief This function will calculate all the meteorological variables based on input
/// @author Edited by XZ Luo
/// @date May 19, 2015


# include "beps.h"


double cal_es(double Ta) {
    return 0.61078 * exp(17.3 * Ta / (237.3 + Ta));
}

double cal_ea(double Ta, double RH) {
    return cal_es(Ta) * RH / 100;
}

double cal_lambda(double Ta) {
    return (2.501 - 0.00237 * Ta) * 1000000;
}

double cal_slope(double Ta) {
    return 2503.0 / pow((Ta + 237.3), 2) * exp(17.27 * Ta / (Ta + 237.3));
}

double ea2q(double ea) {
    return 0.622 * ea / (101.35 - 0.378 * ea);
}


/// @brief Function to calculate meteorological variables based on input
/// @details default input is temperature (C) and relative humidity (0-100)
///          output is an array, named as meteo_pack_output []
/// @details [input] meteo_pack_output [1]= air_density kg/m3 \n
///                  meteo_pack_output [2]= specific heat of air J/kg/C \n
///                  meteo_pack_output [3]= VPD kPa \n
///                  meteo_pack_output [4]= slope of vapor pressure to temperature kPa/C \n
///                  meteo_pack_output [5]= psychrometer constant kPa/C \n
///                  meteo_pack_output [6]= saturate water vapor potential kPa \n
///                  meteo_pack_output [7]= actual water vapor potential kPa \n
///                  meteo_pack_output [8]= specific humidity g/g
/// @param Ta               temperature
/// @param RH                 relative humidity
/// @param meteo_pack_output  meteorological variables array
/// @return void
void meteo_pack(double Ta, double RH, double* meteo_pack_output)
{
    double density_air, cp_air, vpd, slope, gamma, es, ea, q;
    double latent_water;  // latent heat of water J/kg
    density_air = 1.292;
    es = 0.61078 * exp(17.3 * Ta / (237.3 + Ta));
    ea = es * RH / 100;
    vpd = es - ea;
    q = ea2q(ea);
    cp_air = 1004.65 * (1 + 0.84 * q);
    slope = cal_slope(Ta);
    latent_water = (2.501 - 0.00237 * Ta) * 1000000;
    gamma = 0.066;
    // gamma = cp_air * 101.13 / (0.622 * latent_water);

    meteo_pack_output[1]= density_air;
    meteo_pack_output[2]= cp_air;
    meteo_pack_output[3]= vpd;
    meteo_pack_output[4]= slope;
    meteo_pack_output[5]= gamma;
    meteo_pack_output[6]= es;
    meteo_pack_output[7]= ea;
    meteo_pack_output[8]= q;
}
