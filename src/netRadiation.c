/// @file netRadiation.c
/// @brief This module will calculate net radiation at both canopy level and leaf level
/// @author Edited by XZ Luo
/// @date May 23, 2015

#include "beps.h"

double cal_Rl(double emiss, double T) {
    double sb_constant = 5.67 / 100000000;  // stephen-boltzman constant
    return emiss * sb_constant * pow(T + 273.15, 4);
}

/// @brief Function to calculate net radiation at canopy level and leaf level
/// @details [input] global solar radiation,
///                  cosine value for solar zenith angle, albedo of leaves
///                  albedo of snow, percentage of snow cover,
///                  leaf area index overstorey and understorey,
///                  temperature of overstorey, understorey and ground (contain snow?)
///                  temperature of air (C), relative humidity (0-100)
/// @details [output] net radiation for canopy, overstorey, understorey and ground;
///                   net radiation on sunlit, shaded leaves of overstorey and understorey.
/// @param  Rs_global           global short radiation
/// @param  CosZs                     cosine value of solar zenith angle
/// @param  temp_o                    temperature of overstorey
/// @param  temp_u                    temperature of understory
/// @param  temp_g                    temperature of ground
/// @param  lai_o                     leaf area index of overstory, without stem
/// @param  lai_u                     leaf area index of understory, without stem
/// @param  lai_os                    leaf area index of overstory, with stem
/// @param  lai_us                    leaf area index of understory, with stem
/// @param  lai                       sunlit leaves LAI with consideration of stem (PAI)
/// @param  clumping                  clumping index
/// @param  temp_air                  air temperature
/// @param  rh                        relative humidity
/// @param  albedo_snow_v             albedo of snow in this step, visible
/// @param  albedo_snow_n             albedo of snow in this step, near infrared
/// @param  percentArea_snow_o        percentage of snow on overstorey (by area)
/// @param  percentArea_snow_u        percentage of snow on understorey (by area)
/// @param  percent_snow_g            percentage of snow on ground (by mass)
/// @param  albedo_v_o                albedo of overstory, visible, not considering snow, decided by land cover
/// @param  albedo_n_o                albedo of overstory, near infrared
/// @param  albedo_v_u                albedo of understory, visible
/// @param  albedo_n_u                albedo of understory, near infrared
/// @param  albedo_v_g                albedo of ground, visible
/// @param  albedo_n_g                albedo of ground, near infrared
///
/// @param  Rn_o                  net radiation on overstorey
/// @param  Rn_u                  net radiation on understorey
/// @param  Rn_g                  net radiation on ground
/// @param  Rn_Leaf                net radiation at the leaf level, for ET calculation
/// @param  Rns_Leaf           net shortwave radiation at leaf level, for GPP calculation
/// @return void
void netRadiation(double Rs_global, double CosZs, double temp_o, double temp_u, double temp_g,
                  double lai_o, double lai_u, double lai_os, double lai_us,  // LAI of overstorey and understorey, with and without stem
                                                                             //   double lai_o_sunlit, double lai_o_shaded, double lai_u_sunlit, double lai_u_shaded,
                  Leaf lai,
                  double clumping, double temp_air, double rh,
                  double albedo_snow_v, double albedo_snow_n, double percentArea_snow_o, double percentArea_snow_u, double percent_snow_g,
                  double albedo_v_o, double albedo_n_o, double albedo_v_u, double albedo_n_u, double albedo_v_g, double albedo_n_g,
                  double* Rn_o, double* Rn_u, double* Rn_g,
                  Leaf* Rn_Leaf,
                  Leaf* Rns_Leaf
                  //   double* Rn_Leaf.o_sunlit, double* Rn_Leaf.o_shaded, double* Rn_Leaf.u_sunlit, double* Rn_Leaf.u_shaded,
                  //   double* Rns_Leaf.o_sunlit, double* Rns_Leaf.o_shaded, double* Rns_Leaf.u_sunlit, double* neRns_Leaf.u_shaded
) {
    double Rns_o, Rns_u, Rns_g;                                            // net short wave radiation on overstorey, understorey and ground
    double Rns_o_dir, Rns_o_df, Rns_u_dir, Rns_u_df, Rns_g_dir, Rns_g_df;  // direct and diffuse part of net solar radiation
    double Rs_dir, Rs_df;                                                  // direct and diffuse radiation on top of the canopy
    Leaf Rnl_Leaf;
    // double Rnl_Leaf.o_sunlit, Rnl_Leaf.o_shaded, Rnl_Leaf.u_sunlit, Rnl_Leaf.u_shaded;

    double Rnl_o, Rnl_u, Rnl_g;  // net long wave radiation
    double Rs_o_dir, Rs_u_dir, Rs_o_df, Rs_u_df;
    double albedo_o, albedo_u, albedo_g;                                                  // albedo of overstorey, understorey and ground (considering snow)
    double albedo_v_os, albedo_n_os, albedo_v_us, albedo_n_us, albedo_v_gs, albedo_n_gs;  // albedo of three parts in visible and NIR band, (considering snow)

    double e_actual;               // saturated and actual water vapor potential
    double emiss_air, emiss_o, emiss_u, emiss_g;  // emissivity of air, overstorey, understorey and ground
    double Rl_air, Rl_o, Rl_u, Rl_g;              // long wave radiation emitted by different part

    double gap_o_dir, gap_u_dir, gap_o_df, gap_u_df;  // gap fraction of direct and diffuse radiation for overstorey and understorey (diffuse used for diffuse solar radiation and longwave radiation)
    double gap_os_df, gap_us_df;                      // like above, considering stem
    double ratio_cloud;                               // a simple ratio to differentiate diffuse and direct radiation

    // calculate albedo of canopy in this step
    albedo_v_os = albedo_v_o * (1 - percentArea_snow_o) + albedo_snow_v * percentArea_snow_o;  // visible, overstory
    albedo_n_os = albedo_n_o * (1 - percentArea_snow_o) + albedo_snow_n * percentArea_snow_o;  // near infrared
    albedo_v_us = albedo_v_u * (1 - percentArea_snow_u) + albedo_snow_v * percentArea_snow_u;  //        , understory
    albedo_n_us = albedo_n_u * (1 - percentArea_snow_u) + albedo_snow_n * percentArea_snow_u;

    albedo_o = 0.5 * (albedo_v_os + albedo_n_os);
    albedo_u = 0.5 * (albedo_v_us + albedo_n_us);

    // calculate albedo of ground in this step
    albedo_v_gs = albedo_v_g * (1 - percent_snow_g) + albedo_snow_v * percent_snow_g;
    albedo_n_gs = albedo_n_g * (1 - percent_snow_g) + albedo_snow_n * percent_snow_g;
    albedo_g = 0.5 * (albedo_v_gs + albedo_n_gs);

    // separate global solar radiation into direct and diffuse one
    if (CosZs < 0.001)  // solar zenith angle small, all diffuse radiation
        ratio_cloud = 0;
    else
        ratio_cloud = Rs_global / (1367 * CosZs);  // Luo2018, A4

    if (ratio_cloud > 0.8)
        Rs_df = 0.13 * Rs_global;  // Luo2018, A2
    else
        Rs_df = (0.943 + 0.734 * ratio_cloud - 4.9 * pow((ratio_cloud), 2) + 1.796 * pow((ratio_cloud), 3) + 2.058 * pow((ratio_cloud), 4)) * Rs_global;  // Luo2018, A2

    Rs_df = clamp(Rs_df, 0, Rs_global);
    Rs_dir = Rs_global - Rs_df;  // Luo2018, A3

    // fraction at each layer of canopy, direct and diffuse. use Leaf only lai here
    gap_o_dir = exp(-0.5 * clumping * lai_o / CosZs);
    gap_u_dir = exp(-0.5 * clumping * lai_u / CosZs);

    // double gap_os_dir = exp(-0.5 * clumping * lai_os / CosZs);  // considering stem
    // double gap_us_dir = exp(-0.5 * clumping * lai_us / CosZs);
    // indicators to describe leaf distribution angles in canopy. slightly related with LAI
    double cosQ_o = 0.537 + 0.025 * lai_o;  // Luo2018, A10, a representative zenith angle for diffuse radiation transmission
    double cosQ_u = 0.537 + 0.025 * lai_u;

    gap_o_df = exp(-0.5 * clumping * lai_o / cosQ_o);
    gap_u_df = exp(-0.5 * clumping * lai_u / cosQ_u);

    gap_os_df = exp(-0.5 * clumping * lai_os / cosQ_o);  // considering stem
    gap_us_df = exp(-0.5 * clumping * lai_us / cosQ_u);

    e_actual = cal_ea(temp_air, rh);

    emiss_air = 1 - exp(-(pow(e_actual * 10.0, (temp_air + 273.15) / 1200.0)));
    emiss_air = clamp(emiss_air, 0.7, 1.0);

    emiss_o = 0.98;
    emiss_u = 0.98;
    emiss_g = 0.96;

    // net short direct radiation on canopy and ground
    if (Rs_global > zero && CosZs > zero)  // only happens in day time, when sun is out
    {
        // Rns_dir_under = Rs_dir * gap_o_dir * (1 - albedo_u); // from o->u
        Rns_o_dir = Rs_dir * ((1 - albedo_o) - (1 - albedo_u) * gap_o_dir);  // dir into dif_under
        Rns_u_dir = Rs_dir * gap_o_dir * ((1 - albedo_u) - (1 - albedo_g) * gap_u_dir);
        Rns_g_dir = Rs_dir * gap_o_dir * gap_u_dir * (1 - albedo_g);
    } else {
        Rns_o_dir = 0;
        Rns_u_dir = 0;
        Rns_g_dir = 0;
    }

    // net short diffuse radiation on canopy and ground
    if (Rs_global > zero && CosZs > zero)  // only happens in day time, when sun is out
    {
        Rns_o_df = Rs_df * ((1 - albedo_o) - (1 - albedo_u) * gap_o_df) +
                   0.21 * clumping * Rs_dir * (1.1 - 0.1 * lai_o) * exp(-CosZs);  // A8
        Rns_u_df = Rs_df * gap_o_df * ((1 - albedo_u) - (1 - albedo_g) * gap_u_df) +
                   0.21 * clumping * Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_u) * exp(-CosZs);  // A9
        Rns_g_df = Rs_df * gap_o_df * gap_u_df * (1 - albedo_g);
    } else {
        Rns_o_df = 0;
        Rns_u_df = 0;
        Rns_g_df = 0;
    }

    // total net shortwave radiation at canopy level
    Rns_o = Rns_o_dir + Rns_o_df;
    Rns_u = Rns_u_dir + Rns_u_df;
    Rns_g = Rns_g_dir + Rns_g_df;

    // net long wave radiation on canopy and ground
    Rl_air = cal_Rl(emiss_air, temp_air);
    Rl_o = cal_Rl(emiss_o, temp_o);
    Rl_u = cal_Rl(emiss_u, temp_u);
    Rl_g = cal_Rl(emiss_g, temp_g);

    Rnl_o = (emiss_o * (Rl_air + Rl_u * (1 - gap_u_df) + Rl_g * gap_u_df) - 2 * Rl_o) *
                (1 - gap_o_df) +
            emiss_o * (1 - emiss_u) * (1 - gap_u_df) * (Rl_air * gap_o_df + Rl_o * (1 - gap_o_df));

    Rnl_u = (emiss_u * (Rl_air * gap_o_df + Rl_o * (1 - gap_o_df) + Rl_g) - 2 * Rl_u) * (1 - gap_u_df) +
            (1 - emiss_g) * ((Rl_air * gap_o_df + Rl_o * (1 - gap_o_df)) * gap_u_df + Rl_u * (1 - gap_u_df)) +
            emiss_u * (1 - emiss_o) * (Rl_u * (1 - gap_u_df) + Rl_g * gap_u_df) * (1 - gap_o_df);

    Rnl_g = emiss_g * ((Rl_air * gap_o_df + Rl_o * (1 - gap_o_df)) * gap_u_df + Rl_u * (1 - gap_u_df)) -
            Rl_g + (1 - emiss_u) * Rl_g * (1 - gap_u_df);

    // total net radiation for overstorey, understorey and ground.
    *Rn_o = Rns_o + Rnl_o;
    *Rn_u = Rns_u + Rnl_u;
    *Rn_g = Rns_g + Rnl_g;

    // leaf level net radiation updated way
    // reference Chen 2012 clumping index paper
    if (Rs_global > zero && CosZs > zero)  // only happens in day time, when sun is out
    {
        Rs_o_dir = 0.5 * Rs_dir / CosZs;
        Rs_o_dir = min(Rs_o_dir, 0.7 * 1362);
        Rs_u_dir = Rs_o_dir;

        Rs_o_df = (Rs_df - Rs_df * gap_os_df) / lai_os + 0.07 * Rs_dir * (1.1 - 0.1 * lai_os) * exp(-CosZs);
        Rs_u_df = (Rs_df * gap_o_df - Rs_df * gap_o_df * gap_us_df) / lai_us  // pay attention to the type of gap fraction used here
                  + 0.05 * Rs_dir * gap_o_dir * (1.1 - 0.1 * lai_us) * exp(-CosZs);
    } else {
        Rs_o_dir = 0;
        Rs_u_dir = 0;
        Rs_o_df = 0;
        Rs_u_df = 0;
    }

    // overstorey sunlit leaves
    Rns_Leaf->o_sunlit = (Rs_o_dir + Rs_o_df) * (1 - albedo_o);  // diffuse
    Rnl_Leaf.o_sunlit = lai.o_sunlit > 0 ? Rnl_o / lai_os : Rnl_o;

    // overstorey shaded leaf
    Rns_Leaf->o_shaded = Rs_o_df * (1 - albedo_o);  // diffuse
    Rnl_Leaf.o_shaded = lai.o_shaded > 0 ? Rnl_o / lai_os : Rnl_o;

    // understorey sunlit leaf
    Rns_Leaf->u_sunlit = (Rs_u_dir + Rs_u_df) * (1 - albedo_u);
    Rnl_Leaf.u_sunlit = lai.u_sunlit > 0 ? Rnl_u / lai_us : Rnl_u;

    // understorey shaded leaf
    Rns_Leaf->u_shaded = Rs_u_df * (1 - albedo_u);
    Rnl_Leaf.u_shaded = lai.u_shaded > 0 ? Rnl_u / lai_us : Rnl_u;

    Rn_Leaf->o_sunlit = Rns_Leaf->o_sunlit + Rnl_Leaf.o_sunlit;
    Rn_Leaf->o_shaded = Rns_Leaf->o_shaded + Rnl_Leaf.o_shaded;
    Rn_Leaf->u_sunlit = Rns_Leaf->u_sunlit + Rnl_Leaf.u_sunlit;
    Rn_Leaf->u_shaded = Rns_Leaf->u_shaded + Rnl_Leaf.u_shaded;
}

//  leaf level net radiation: original way, use canopy radiation divided by LAI.
//  here calculate net radiation for all leaf component again use the original method.
//  These code could be commented, because not right, I just put it here to validate model.
//    if(lai_o_sunlit>0)
//        Rn_Leaf->o_sunlit = Rns_o_dir/lai_o_sunlit + Rnl_o/lai_os;
//    else
//        Rn_Leaf->o_sunlit = Rns_o_dir + Rnl_o;
//
//    if(lai_o_shaded>0)
//        Rn_Leaf->o_shaded = Rns_o_df/lai_o_shaded + Rnl_o/lai_os;
//    else
//        Rn_Leaf->o_shaded = Rns_o_df + Rnl_o;
//
//
//    if(lai_u_sunlit>0)
//        Rn_Leaf->u_sunlit = Rns_u_dir/lai_u_sunlit + Rnl_u/lai_us;
//    else
//        Rn_Leaf->u_sunlit = Rns_u_dir + Rnl_u;
//
//
//    if(lai_u_shaded>0)
//        Rn_Leaf->u_shaded = Rns_u_df/lai_u_shaded + Rnl_u/lai_us;
//    else
//        Rn_Leaf->u_shaded =Rns_u_df + Rnl_u;
