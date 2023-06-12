/// @file inter_prg.c
/// @brief the inter-program between main program and modules
/// @date Last update: July, 2015

#include "beps.h"
#include "soil.h"

/// @brief the inter-module function between main program and modules
/// @param  jday       day of year
/// @param  rstep      hour of day
/// @param  lai        leaf area index
/// @param  clumping   clumping index
/// @param  parameter  parameter array according to land cover types
/// @param  meteo      meteorological data
/// @param  CosZs      cosine of solar zenith angle
/// @param  var_o      temporary variables array of last time step
/// @param  var_n      temporary variables array of this time step
/// @param  soilp      soil coefficients according to land cover types and soil textures
/// @param  mid_res    results struct
/// @return void
void inter_prg_c(int jday, int rstep, double lai, double clumping, double parameter[], struct climatedata* meteo,
                 double CosZs, double var_o[], double var_n[], struct Soil* soilp, struct results* mid_res, struct OutputET* mid_ET) {
    /*****  define parameters and arrays  *****/
    int num, kkk;
    int landcover;
    double lai_o, lai_u;
    double stem_o, stem_u;

    double d_soil[layer + 1];
    double Zsp;  // the depth of snow on the surface
    double Zp;   // depth of pounded water on the surface
    // double Zp1 = 0, makeZp2 = 0;
    double height_wind_sp;  // height of the Va measured for calculation of L

    double Qhc_o[MAX_Loop] = {0}, Qhc_u[MAX_Loop] = {0}, Qhg[MAX_Loop] = {0};  // The sensible heat flux from canopy and ground
    double G[layer + 2][MAX_Loop] = {0};                                       // the heat flux into the canopy of over story --in W/m^2

    double Wcl_o[MAX_Loop] = {0}, Wcs_o[MAX_Loop] = {0};  // the masses od rain and snow on the canopy
    double Xcl_o[MAX_Loop] = {0}, Xcs_o[MAX_Loop] = {0};  // the fraction of canopy covered by liquid water and snow
    double Wcl_u[MAX_Loop] = {0}, Wcs_u[MAX_Loop] = {0};  // the masses of rain and snow on the canopy
    double Xcl_u[MAX_Loop] = {0}, Xcs_u[MAX_Loop] = {0};  // the fraction of canopy covered by liquid water and snow

    double r_rain_g[MAX_Loop] = {0};                                // the rainfall rate, on ground surface  m/s
    double rho_snow[MAX_Loop] = {0};                                // density of snow
    double alpha_v_sw[MAX_Loop] = {0}, alpha_n_sw[MAX_Loop] = {0};  // albedo of snow
    double Wg_snow[MAX_Loop] = {0};                                 // the amount of snow on the ground
    double Xg_snow[MAX_Loop] = {0};                                 // the fraction of the ground surface covered by snow
    double Ac_snow_u[MAX_Loop] = {0};                               // the areas of canopy covered in snow, o--overstory and  u-- understory
    double Ac_snow_o[MAX_Loop] = {0};

    double Ts0[MAX_Loop] = {0}, Tsn0[MAX_Loop] = {0}, Tsm0[MAX_Loop] = {0}, Tsn1[MAX_Loop] = {0}, Tsn2[MAX_Loop] = {0};  // surface temperature
    double Tc_u[MAX_Loop] = {0};                                                                                         // the effective canopy temperature in K
    double Tm[layer + 2][MAX_Loop] = {0};                                                                                // Tb[layer+2][MAX_Loop],soil temperature at the bottom and the middle of each layer*/

    double lambda[layer + 2] = {0};        // thermal conductivity of each soil layer;*/
    double Cs[layer + 2][MAX_Loop] = {0};  // the soil volumetric heat capacity of each soil layer, j+kkk/m^3/K*/

    double temp_air;                 // air temperature */
    double precip, rh_air, wind_sp;  // precipitation in meters, relative humidity (%), wind speed in m/s */
    double temp_grd;                 // ground temperature */

    double Eil_o[MAX_Loop] = {0}, EiS_o[MAX_Loop] = {0};      // the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snow
    double Eil_u[MAX_Loop] = {0}, EiS_u[MAX_Loop] = {0};      // the evaporation rate of intercepted water of overstory--in kg/m^2/s; l-- water; S-snow of intercepted water--in kg/m^2/s
    double Trans_o[MAX_Loop] = {0}, Trans_u[MAX_Loop] = {0};  // transpiration from overstory and understory canopies
    double Evap_soil[MAX_Loop] = {0};                         // evaporation from soil
    double Evap_SW[MAX_Loop] = {0};                           // evaporation from water pond
    double Evap_SS[MAX_Loop] = {0};                           // evaporation from snow pack

    double lambda_snow[MAX_Loop] = {0};  // the effective thermal conductivity of snow --in m^2/s
    double e_a10;                        // the vapour partial pressure of water --in kPa(1mb=100Pa=0.1kpa)

    double Ks;  // KsCal,KsMea[MAX_Loop] instantaneous total short wave radiation (Global radiation)
    double alpha_sat, alpha_dry;
    double alpha_v_o, alpha_v_u;  // visible albedo of overstory,  o--overstory, u--understory;
    double alpha_n_o, alpha_n_u;  // near-infrared albedo of overstory,o--overstory, u--understory;
    double alpha_g;               // the all-wave ground surface albedo
    double alpha_v_g, alpha_n_g;  // the ground surface albedo for visible range and for near-infrared respectively
    double Cp_ca;                 // specific heat of moist air above the canopy  in j+kkk/k/G
    double ra_o, ra_u, ra_g;      // the aerodynamic resistance of overstory, understory and ground surface   in s/m

    double q_ca;  // the actual canopy stomatal resistance  --in s/m
    double Rn_o, Rn_u, Rn_g;

    // double ip = 0;  // the cumulative infiltration at the time of ponding   --in m/s
    // double Inf = 0;
    // double zr = 0.8;
    double Cpd = 1004.65;
    double rho_w = 1025.0;  // density of water

    Leaf Cc_new;  // CO2 concentration in the chloroplast
    Leaf Cs_old;  // CO2 concentration on the surfaces of leaves
    Leaf Cs_new;
    Leaf Ci_old;  // intercellular CO2 concentration pn the leaf
    Leaf Ci_new;
    Leaf Tc_old;  // the effective canopy temperature in K
    Leaf Tc_new;
    Leaf Gs_old;  // stomatal conductance of the big leaf for water
    Leaf Gs_new;
    Leaf Gc;   // the total conductance for CO2 from the intercellular space of the leaves to the reference height above the canopy
    Leaf Gw;   // the total conductance for water from the intercellular space of the leaves to the reference height above the canopy
    Leaf Gww;  // the total conductance for water from the surface of the leaves to the reference height above the canopy
    Leaf Gh;   // total conductance for heat transfer from the leaf surface to the reference height above the canopy
    Leaf Ac;   // net photosynthesis rate

    Leaf Rn_Leaf;  // net radiation of leaves
    Leaf Rns_Leaf;          // Rabs, absorbed solar radiation
    Leaf leleaf;     // leaf latent heat flux (mol m-2 s-1)
    Leaf GPP, LAI;
    Leaf PAI;  // PAI = LAI + SAI; // double LAI.o_sunlit, LAI.o_shaded, LAI.u_sunlit, LAI.u_shaded;

    double f_soilwater;  // an empirical parameter describing the relative availability of soil water for plants
    double psychrometer = 0.066;
    double Tco, Tcu, slope;
    double H_o_sunlit, H_o_shaded;  // sensible heat flux from leaves

    double VPS_air;
    double GH_o, Ga_o, Gb_o, Ga_u, Gb_u;
    double canopyh_o, canopyh_u;

    double VPD_air;  // Vapor pressure deficit term

    double mass_water_g;
    double percentArea_snow_o, percentArea_snow_u;
    double Gheat_g;

    double b_h2o;  // the intercept term in BWB model (mol H2O m-2 s-1)
    double m_h2o;  // the slope in BWB model

    // leaf latent heat flux (mol/m2/s)
    // double leleaf.o_sunlit, leleaf.o_shaded, leleaf.u_sunlit, leleaf.u_shaded;

    // parameters for Vcmax-Nitrogen calculations
    double Kn = 0.3;  // 0.713/2.4
    double G_theta = 0.5;
    double K, Vcmax0, Vcmax_sunlit, Vcmax_shaded, expr1, expr2, expr3;
    // double slope_Vcmax_N, leaf_N, Jmax_sunlit, Jmax_shaded;

    alpha_sat = parameter[24];  // albedo of saturated/dry soil for module rainfall 1
    alpha_dry = parameter[25];  // the albedo of dry soil
    canopyh_o = parameter[29];  // to be used for module aerodynamic_conductance
    canopyh_u = parameter[30];
    height_wind_sp = parameter[31];  // the_height_to_measure_wind_speed, for module aerodynamic_conductance
    m_h2o = parameter[33];           // to be used for module photosynthesis
    b_h2o = parameter[34];

    /*****  Vcmax-Nitrogen calculations，by G.Mo，Apr. 2011  *****/
    if (CosZs > 0)  // day time
    {
        K = G_theta * clumping / CosZs;  // G_theta = 0.5 assuming a spherical leaf angle distribution
        Vcmax0 = parameter[36];
        expr1 = 1 - exp(-K * lai);
        expr2 = 1 - exp(-lai * (Kn + K));
        expr3 = 1 - exp(-Kn * lai);

        // Formulas based on Chen et al., 2012, GBC
        if (expr1 > 0)
            Vcmax_sunlit = Vcmax0 * parameter[47] * parameter[46] * K * expr2 / (Kn + K) / expr1;
        else
            Vcmax_sunlit = Vcmax0;

        if (K > 0 && lai > expr1 / K)
            Vcmax_shaded = Vcmax0 * parameter[47] * parameter[46] * (expr3 / Kn - expr2 / (Kn + K)) / (lai - expr1 / K);
        else
            Vcmax_shaded = Vcmax0;
    }

    /*****  LAI calculation module, by B. Chen  *****/
    lai_o = lai;
    if (lai < 0.1) lai_o = 0.1;
    landcover = (int)parameter[4];

    if (landcover == 25 || landcover == 40)
        lai_u = 0.01;
    else
        lai_u = 1.18 * exp(-0.99 * lai_o);

    if (lai_u > lai_o) lai_u = 0.01;

    stem_o = parameter[8] * 0.2;  // parameter[8]->LAI max overstory
    stem_u = parameter[9] * 0.2;  // parameter[9]->LAI max understory

    // lai_calc module
    // separate lai into sunlit and shaded portions
    lai2(clumping, CosZs, stem_o, stem_u, lai_o, lai_u, &LAI, &PAI);
    /*****  Initialization of this time step  *****/

    Ks = meteo->Srad;
    rh_air = meteo->rh;
    wind_sp = meteo->wind;
    precip = meteo->rain / step;  // precipitation in meters
    temp_air = meteo->temp;

    if (Ks <= 0) {
        alpha_v_o = 0;
        alpha_n_o = 0;
        alpha_v_u = 0;
        alpha_n_u = 0;
    } else {
        alpha_v_o = parameter[22];
        alpha_n_o = parameter[23];
        alpha_v_u = parameter[22];
        alpha_n_u = parameter[23];
    }

    // Ground surface temperature
    Ts0[0] = clamp(var_o[3], temp_air - 2.0, temp_air + 2.0);
    Tsn0[0] = clamp(var_o[4], temp_air - 2.0, temp_air + 2.0);
    Tsm0[0] = clamp(var_o[5], temp_air - 2.0, temp_air + 2.0);
    Tsn1[0] = clamp(var_o[6], temp_air - 2.0, temp_air + 2.0);
    Tsn2[0] = clamp(var_o[7], temp_air - 2.0, temp_air + 2.0);

    Qhc_o[0] = var_o[11];
    Wcl_o[0] = var_o[15];
    Wcs_o[0] = var_o[16]; /* the mass of intercepted liquid water and snow, overstory */

    // the evaporation rate of rain and snow--in kg/m^2/s, understory
    Wcl_u[0] = var_o[18];
    Wcs_u[0] = var_o[19]; /* the mass of intercepted liquid water and snow, overstory */

    Wg_snow[0] = var_o[20]; /* thr fraction of ground surface covered in snow and snow mass */

    Zsp = soilp->Zsp;
    Zp = soilp->Zp;

    if (Zp < 0.001) Zp = 0;
    /*****  Vcmax Jmax module by L. He  *****/
    // slope_Vcmax_N = parameter[47];
    // leaf_N = parameter[46];

    // Vcmax_Jmax(lai_o, clumping, Vcmax0,slope_Vcmax_N, leaf_N, CosZs, &Vcmax_sunlit, &Vcmax_shaded, &Jmax_sunlit, &Jmax_shaded);

    // temperatures of overstorey and understorey canopies
    init_leaf_dbl(&Tc_old, temp_air - 0.5);

    /*****  Ten time intervals in a hourly time step->6min or 360s per loop  ******/
    for (kkk = 1; kkk <= kloop; kkk++) {
        /*****  Snow pack stage 1 by X. Luo  *****/
        snowpack_stage1(temp_air, precip, Wcs_o[kkk - 1], Wcs_u[kkk - 1], Wg_snow[kkk - 1],
                        &Wcs_o[kkk], &Wcs_u[kkk], &Wg_snow[kkk], lai_o, lai_u, clumping,
                        &Ac_snow_o[kkk], &Ac_snow_u[kkk], &Xcs_o[kkk], &Xcs_u[kkk], &Xg_snow[kkk],
                        &rho_snow[kkk], &Zsp, &alpha_v_sw[kkk], &alpha_n_sw[kkk]);

        /*****  Rain fall stage 1 by X. Luo  *****/
        rainfall_stage1(temp_air, precip, Wcl_o[kkk - 1], Wcl_u[kkk - 1],
                        lai_o, lai_u, clumping, &Wcl_o[kkk], &Wcl_u[kkk], &Xcl_o[kkk], &Xcl_u[kkk], &r_rain_g[kkk]);

        // Old version
        // if(thetam[0][kkk-1]<soilp->theta_vwp[1]*0.5) alpha_g = alpha_dry;
        // else alpha_g = (thetam[0][kkk-1]-soilp->theta_vwp[1]*0.5)/(soilp->fei[1]-soilp->theta_vwp[1]*0.5) * (alpha_sat - alpha_dry) + alpha_dry;
        if (soilp->thetam_prev[1] < soilp->theta_vwp[1] * 0.5)
            alpha_g = alpha_dry;
        else
            alpha_g = (soilp->thetam_prev[1] - soilp->theta_vwp[1] * 0.5) / (soilp->fei[1] - soilp->theta_vwp[1] * 0.5) * (alpha_sat - alpha_dry) + alpha_dry;

        alpha_v_g = 2.0 / 3.0 * alpha_g;
        alpha_n_g = 4.0 / 3.0 * alpha_g;

        /*****  Soil water factor module by L. He  *****/
        soil_water_factor_v2(soilp);
        f_soilwater = soilp->f_soilwater;

        if (f_soilwater > 1.0) f_soilwater = 1.0;  // to be used for module photosynthesis

        GH_o = Qhc_o[kkk - 1];  // to be used as the init. for module aerodynamic_conductance

        VPS_air = 0.61078 * exp(17.3 * temp_air / (237.3 + temp_air));  // to estimate saturated water vapor pressure in kpa
        e_a10 = VPS_air * rh_air / 100;                                 // to be used for module photosynthesis
        VPD_air = VPS_air - e_a10;                                      // water vapor deficit at the reference height

        q_ca = 0.622 * e_a10 / (101.35 - 0.378 * e_a10);  // in g/g, unitless
        Cp_ca = Cpd * (1 + 0.84 * q_ca);

        slope = 2503.0 / pow((temp_air + 237.3), 2) * exp(17.27 * temp_air / (temp_air + 237.3));

        init_leaf_dbl(&Ci_old, 0.7 * CO2_air);          // over- and under-store
        init_leaf_dbl2(&Gs_old, 1 / 200.0, 1 / 300.0);  // over- and under-store

        percentArea_snow_o = Ac_snow_o[kkk] / lai_o / 2;
        percentArea_snow_u = Ac_snow_u[kkk] / lai_u / 2;

        temp_grd = temp_air;  // ground temperature substituted by air temperature

        num = 0;
        while (1)  // iteration for BWB equation until results converge
        {
            num = num + 1;

            /***** Aerodynamic_conductance module by G.Mo  *****/
            aerodynamic_conductance(canopyh_o, canopyh_u, height_wind_sp, clumping, temp_air, wind_sp, GH_o,
                                    lai_o + stem_o, lai_u + stem_u, &ra_o, &ra_u, &ra_g, &Ga_o, &Gb_o, &Ga_u, &Gb_u);

            Gh.o_sunlit = 1.0 / (1.0 / Ga_o + 0.5 / Gb_o);  // heat conductance of sunlit leaves of overstorey
            Gh.o_shaded = 1.0 / (1.0 / Ga_o + 0.5 / Gb_o);  // heat conductance of shaded leaves of overstorey
            Gh.u_sunlit = 1.0 / (1.0 / Ga_u + 0.5 / Gb_u);  // heat conductance of sunlit leaves of understorey
            Gh.u_shaded = 1.0 / (1.0 / Ga_u + 0.5 / Gb_u);  // heat conductance of shaded leaves of understorey

            Gww.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 100);  // conductance for intercepted water of sunlit leaves of overstorey
            Gww.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 100);  // conductance for intercepted water of shaded leaves of overstorey
            Gww.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 100);  // conductance for intercepted water of sunlit leaves of understorey
            Gww.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 100);  // conductance for intercepted water of shaded leaves of understorey

            // temperatures of overstorey and understorey canopies
            Tco = (Tc_old.o_sunlit * PAI.o_sunlit + Tc_old.o_shaded * PAI.o_shaded) / (PAI.o_sunlit + PAI.o_shaded);
            Tcu = (Tc_old.u_sunlit * PAI.u_sunlit + Tc_old.u_shaded * PAI.u_shaded) / (PAI.u_sunlit + PAI.u_shaded);

            /*****  Net radiation at canopy and leaf level module by X.Luo  *****/
            netRadiation(Ks, CosZs, Tco, Tcu, temp_grd, lai_o, lai_u, lai_o + stem_o, lai_u + stem_u,
                         PAI,
                         clumping, temp_air, rh_air, alpha_v_sw[kkk], alpha_n_sw[kkk],
                         percentArea_snow_o, percentArea_snow_u,
                         Xg_snow[kkk], alpha_v_o, alpha_n_o, alpha_v_u, alpha_n_u, alpha_v_g, alpha_n_g,
                         &Rn_o, &Rn_u, &Rn_g,
                         &Rn_Leaf, &Rns_Leaf);

            /*****  Photosynthesis module by B. Chen  *****/
            // Four components: overstory sunlit and shaded, understory sunlit and shaded
            Gw.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_old.o_sunlit);  // conductance of sunlit leaves of overstorey for water
            Gw.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_old.o_shaded);  // conductance of shaded leaves of overstorey for water
            Gw.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_old.u_sunlit);  // conductance of sunlit leaves of understorey for water
            Gw.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_old.u_shaded);  // conductance of shaded leaves of understorey for water

            leleaf.o_sunlit = Gw.o_sunlit * (VPD_air + slope * (Tc_old.o_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer;
            leleaf.o_shaded = Gw.o_shaded * (VPD_air + slope * (Tc_old.o_shaded - temp_air)) * rho_a * Cp_ca / psychrometer;
            leleaf.u_sunlit = Gw.u_sunlit * (VPD_air + slope * (Tc_old.u_sunlit - temp_air)) * rho_a * Cp_ca / psychrometer;
            leleaf.u_shaded = Gw.u_shaded * (VPD_air + slope * (Tc_old.u_shaded - temp_air)) * rho_a * Cp_ca / psychrometer;

            if (CosZs > 0) {
                photosynthesis(Tc_old.o_sunlit, Rns_Leaf.o_sunlit, e_a10, Gb_o, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o, Ci_old.o_sunlit,
                               temp_air, leleaf.o_sunlit, &Gs_new.o_sunlit, &Ac.o_sunlit, &Ci_new.o_sunlit);
                photosynthesis(Tc_old.o_shaded, Rns_Leaf.o_shaded, e_a10, Gb_o, Vcmax_shaded, f_soilwater, b_h2o, m_h2o, Ci_old.o_shaded,
                               temp_air, leleaf.o_shaded, &Gs_new.o_shaded, &Ac.o_shaded, &Ci_new.o_shaded);
                photosynthesis(Tc_old.u_sunlit, Rns_Leaf.u_sunlit, e_a10, Gb_u, Vcmax_sunlit, f_soilwater, b_h2o, m_h2o, Ci_old.u_sunlit,
                               temp_air, leleaf.u_sunlit, &Gs_new.u_sunlit, &Ac.u_sunlit, &Ci_new.u_sunlit);
                photosynthesis(Tc_old.u_shaded, Rns_Leaf.u_shaded, e_a10, Gb_u, Vcmax_shaded, f_soilwater, b_h2o, m_h2o, Ci_old.u_shaded,
                               temp_air, leleaf.u_shaded, &Gs_new.u_shaded, &Ac.u_shaded, &Ci_new.u_shaded);
            } else {
                init_leaf_dbl(&Gs_new, 0.0001);
                init_leaf_dbl(&Ac, 0.0);
                init_leaf_dbl(&Cs_new, CO2_air);
                init_leaf_dbl(&Ci_new, CO2_air * 0.7);
                init_leaf_dbl(&Cc_new, CO2_air * 0.7 * 0.8);  // not used
            }

            // init_leaf_struct(&Ci_old, Ci_new);
            // init_leaf_struct(&Cs_old, Cs_new);
            // init_leaf_struct(&Gs_old, Gs_new);  // m/s

            Ci_old.o_sunlit = Ci_new.o_sunlit;
            Cs_old.o_sunlit = Cs_new.o_sunlit;
            Gs_old.o_sunlit = Gs_new.o_sunlit;  // m/s

            Ci_old.o_shaded = Ci_new.o_shaded;
            Cs_old.o_shaded = Cs_new.o_shaded;
            Gs_old.o_shaded = Gs_new.o_shaded;  // m/s

            // note error at here, kongdd, v20220705
            Ci_old.u_sunlit = Ci_new.o_sunlit;
            Cs_old.u_sunlit = Cs_new.u_sunlit;
            Gs_old.u_sunlit = Gs_new.o_sunlit;  // m/s

            Ci_old.u_shaded = Ci_new.u_shaded;
            Cs_old.u_shaded = Cs_new.u_shaded;
            Gs_old.u_shaded = Gs_new.u_shaded;  // m/s

            Gw.o_sunlit = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_sunlit);  // conductance for water
            Gw.o_shaded = 1.0 / (1.0 / Ga_o + 1.0 / Gb_o + 1.0 / Gs_new.o_shaded);
            Gw.u_sunlit = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_sunlit);
            Gw.u_shaded = 1.0 / (1.0 / Ga_u + 1.0 / Gb_u + 1.0 / Gs_new.u_shaded);

            Gc.o_sunlit = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_sunlit);  // conductance for CO2
            Gc.o_shaded = 1.0 / (1.0 / Ga_o + 1.4 / Gb_o + 1.6 / Gs_new.o_shaded);
            Gc.u_sunlit = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_sunlit);
            Gc.u_shaded = 1.0 / (1.0 / Ga_u + 1.4 / Gb_u + 1.6 / Gs_new.u_shaded);

            /***** Leaf temperatures module by L. He  *****/
            Leaf_Temperatures(temp_air, slope, psychrometer, VPD_air, Cp_ca,
                              Gw, Gww, Gh,
                              Xcs_o[kkk], Xcl_o[kkk], Xcs_u[kkk], Xcl_u[kkk],
                              Rn_Leaf,
                              &Tc_new);

            H_o_sunlit = (Tc_new.o_sunlit - temp_air) * rho_a * Cp_ca * Gh.o_sunlit;
            H_o_shaded = (Tc_new.o_shaded - temp_air) * rho_a * Cp_ca * Gh.o_shaded;
            GH_o = H_o_sunlit * PAI.o_sunlit + H_o_shaded * PAI.o_shaded;  // for next num aerodynamic conductance calculation

            if (fabs(Tc_new.o_sunlit - Tc_old.o_sunlit) < 0.02 && fabs(Tc_new.o_shaded - Tc_old.o_shaded) < 0.02 &&
                fabs(Tc_new.u_sunlit - Tc_old.u_sunlit) < 0.02 && fabs(Tc_new.u_shaded - Tc_old.u_shaded) < 0.02)

                break;          // break the iteration if results converge
            else if (num > 22)  // if the iteration does not converge
            {
                init_leaf_dbl(&Tc_old, temp_air);
                break;
            } else {
                init_leaf_struct(&Tc_old, Tc_new);
            }

        }  // end of while

        GPP.o_sunlit = Ac.o_sunlit * LAI.o_sunlit;
        GPP.o_shaded = Ac.o_shaded * LAI.o_shaded;
        GPP.u_sunlit = Ac.u_sunlit * LAI.u_sunlit;
        GPP.u_shaded = Ac.u_shaded * LAI.u_shaded;

        /*****  Transpiration module by X. Luo  *****/
        // transpiration(Tc_new.o_sunlit, Tc_new.o_shaded, Tc_new.u_sunlit, Tc_new.u_shaded, temp_air, rh_air,
        //               Gw.o_sunlit, Gw.o_shaded, Gw.u_sunlit, Gw.u_shaded, LAI.o_sunlit, LAI.o_shaded, LAI.u_sunlit, LAI.u_shaded,
        //               &Trans_o[kkk], &Trans_u[kkk]);
        transpiration(
            Tc_new, temp_air, rh_air,
            Gw, LAI,
            &Trans_o[kkk], &Trans_u[kkk]);

        /*****  Evaporation and sublimation from canopy by X. Luo  *****/
        evaporation_canopy(
            Tc_new, temp_air, rh_air,
            Gww, PAI,
            Xcl_o[kkk], Xcl_u[kkk], Xcs_o[kkk], Xcs_u[kkk],
            &Eil_o[kkk], &Eil_u[kkk], &EiS_o[kkk], &EiS_u[kkk]);

        /*****  Rainfall stage 2 by X. Luo  *****/
        rainfall_stage2(Eil_o[kkk], Eil_u[kkk], &Wcl_o[kkk], &Wcl_u[kkk]);

        /*****  Snow pack stage 2 by X. Luo  *****/
        snowpack_stage2(EiS_o[kkk], EiS_u[kkk], &Wcs_o[kkk], &Wcs_u[kkk]);

        /*****  Evaporation from soil module by X. Luo  *****/
        Gheat_g = 1 / ra_g;
        mass_water_g = rho_w * Zp;

        evaporation_soil(temp_grd, Ts0[kkk - 1], rh_air, Rn_g, Gheat_g, &Xg_snow[kkk],
                         &Zp, &Zsp, &mass_water_g, &Wg_snow[kkk], rho_snow[kkk], soilp->thetam_prev[0], soilp->fei[0],
                         &Evap_soil[kkk], &Evap_SW[kkk], &Evap_SS[kkk]);

        /*****  Soil Thermal Conductivity module by L. He  *****/
        UpdateSoilThermalConductivity(soilp);
        Update_Cs(soilp);

        /*****  Surface temperature by X. Luo  *****/
        Cs[0][kkk] = soilp->Cs[0];  // added
        Cs[1][kkk] = soilp->Cs[0];
        Tc_u[kkk] = Tcu;  // added
        lambda[1] = soilp->lambda[0];
        d_soil[1] = soilp->d_soil[0];
        Tm[1][kkk - 1] = soilp->temp_soil_p[1];  // first place is temp_soil_p[0]?
        Tm[0][kkk - 1] = soilp->temp_soil_p[0];
        G[1][kkk] = soilp->G[0];

        surface_temperature(temp_air, rh_air, Zsp, Zp,
                            Cs[1][kkk], Cs[0][kkk], Gheat_g, d_soil[1], rho_snow[kkk], Tc_u[kkk],
                            Rn_g, Evap_soil[kkk], Evap_SW[kkk], Evap_SS[kkk],
                            lambda[1], Xg_snow[kkk], G[1][kkk],
                            Ts0[kkk - 1], Tm[1][kkk - 1], Tm[0][kkk - 1], Tsn0[kkk - 1],
                            Tsm0[kkk - 1], Tsn1[kkk - 1], Tsn2[kkk - 1],
                            &Ts0[kkk], &Tm[0][kkk], &Tsn0[kkk],
                            &Tsm0[kkk], &Tsn1[kkk], &Tsn2[kkk],
                            &G[0][kkk]);
        soilp->temp_soil_c[0] = Tm[0][kkk];

        /*****  Snow Pack Stage 3 module by X. Luo  *****/
        snowpack_stage3(temp_air, Tsn0[kkk], Tsn0[kkk - 1], rho_snow[kkk], &Zsp, &Zp, &Wg_snow[kkk]);

        /*****  Sensible heat flux module by X. Luo  *****/
        sensible_heat(Tc_new, Ts0[kkk], temp_air, rh_air,
                      Gh, Gheat_g, PAI, &Qhc_o[kkk], &Qhc_u[kkk], &Qhg[kkk]);

        // printf("===========================\n");
        // printf("kkk = %d\n", kkk+1);
        // printf("Ts0 = %f, Tm = %f, Tsno = %f, Tsm0 = %f, Tsn1 = %f, Tsn2 = %f, G = %f\n",
        //        Ts0[kkk], Tm[0][kkk], Tsn0[kkk],
        //        Tsm0[kkk], Tsn1[kkk], Tsn2[kkk],
        //        G[0][kkk]);
        // printf("Ta = %f, Gg = %f, QHs = %f, %f, %f\n",
        //        temp_air, Gheat_g,
        //        Qhc_o[kkk], Qhc_u[kkk], Qhg[kkk]);

        /*****  Soil water module by L. He  *****/
        // in-process value check
        // if(jday==75 && rstep==8 && kkk==5)
        // if(jday==75 && rstep==8 )
        // printf("%d, %d, %f, %f\n", jday, rstep, soilp->thetam_prev[0], soilp->f_soilwater);
        soilp->Zsp = Zsp;
        soilp->G[0] = G[0][kkk];

        UpdateHeatFlux(soilp, Xg_snow[kkk], lambda_snow[kkk], Tsn0[kkk], temp_air, kstep);
        Soil_Water_Uptake(soilp, Trans_o[kkk], Trans_u[kkk], Evap_soil[kkk]);

        soilp->r_rain_g = r_rain_g[kkk];
        soilp->Zp = Zp;

        UpdateSoilMoisture(soilp, kstep);
        Zp = soilp->Zp;
    }  // The end of kkk loop

    kkk = kloop;  // the last step

    Tsn1[kkk] = clamp(Tsn1[kkk], -40., 40.);
    Tsn2[kkk] = clamp(Tsn2[kkk], -40., 40.);
    // if (Tsn1[kkk] > 40) Tsn2[kkk] = 40;
    // if (Tsn1[kkk] < -40) Tsn2[kkk] = -40;
    // if (Tsn2[kkk] > 40) Tsn2[kkk] = 40;
    // if (Tsn2[kkk] < -40) Tsn2[kkk] = -40;

    mid_ET->Trans_o   = Trans_o[kkk];
    mid_ET->Trans_u   = Trans_u[kkk];
    mid_ET->Eil_o     = Eil_o[kkk];
    mid_ET->Eil_u     = Eil_u[kkk];
    mid_ET->EiS_o     = EiS_o[kkk];
    mid_ET->EiS_u     = EiS_u[kkk];
    mid_ET->Evap_soil = Evap_soil[kkk];
    mid_ET->Evap_SW   = Evap_SW[kkk];
    mid_ET->Evap_SS   = Evap_SS[kkk];
    mid_ET->Qhc_o     = Qhc_o[kkk];
    mid_ET->Qhc_u     = Qhc_u[kkk];
    mid_ET->Qhg       = Qhg[kkk];

    update_ET(mid_ET, mid_res, temp_air);

    var_n[3] = Ts0[kkk];   // To: The temperature of ground surface
    var_n[4] = Tsn0[kkk];  // To: The temperature of ground surface
    var_n[5] = Tsm0[kkk];
    var_n[6] = Tsn1[kkk];  // To: The temperature of ground surface
    var_n[7] = Tsn2[kkk];  // To: The temperature of ground surface

    var_n[11] = Qhc_o[kkk];
    var_n[15] = Wcl_o[kkk];
    var_n[16] = Wcs_o[kkk];  // the mass of intercepted liquid water and snow, overstory
    var_n[18] = Wcl_u[kkk];
    var_n[19] = Wcs_u[kkk];    // the mass of intercepted liquid water and snow, overstory
    var_n[20] = Wg_snow[kkk];  // the fraction of ground surface covered by snow and snow mass

    mid_res->Net_Rad = Rn_o + Rn_u + Rn_g;
    mid_res->gpp_o_sunlit = GPP.o_sunlit;  // umol C/m2/s
    mid_res->gpp_u_sunlit = GPP.u_sunlit;
    mid_res->gpp_o_shaded = GPP.o_shaded;
    mid_res->gpp_u_shaded = GPP.u_shaded;
    // total GPP -> gC/m2/step
    mid_res->GPP = (GPP.o_sunlit + GPP.o_shaded + GPP.u_sunlit + GPP.u_shaded) * 12 * step * 0.000001;
    return;
}

void update_ET(struct OutputET* x, struct results* mid_res, double Ta) {
    double Lv_liquid = (2.501 - 0.00237 * Ta) * 1000000;  // The latent heat of water vaporization in j/kg
    double Lv_solid = 2.83 * 1000000;  // the latent heat of evaporation from solid (snow/ice) at air temperature=Ta, in j+kkk/kg

    x->LH = Lv_liquid * (x->Trans_o + x->Eil_o + x->Trans_u + x->Eil_u + x->Evap_soil + x->Evap_SW) +
            Lv_solid * (x->EiS_o + x->EiS_u + x->Evap_SS);
    //  LH:  total latent heat flux
    x->SH = x->Qhc_o + x->Qhc_u + x->Qhg;  // SH: total sensible heat flux

    x->Trans = (x->Trans_o + x->Trans_u) * step;                                                            // total transpiration  mm/step
    x->Evap = (x->Eil_o + x->Eil_u + x->Evap_soil + x->Evap_SW + x->EiS_o + x->EiS_u + x->Evap_SS) * step;  // total evaporation -> mm/step

    mid_res->LH = x->LH;
    mid_res->SH    = x->SH;
    mid_res->Trans = x->Trans;
    mid_res->Evap = x->Evap;
}

void inter_prg(int jday, int rstep, double lai, double clumping, double parameter[], struct climatedata* meteo,
               double CosZs, double var_o[], double var_n[], struct Soil* p_soil, struct results* mid_res) {
    struct OutputET* mid_ET;
    inter_prg_c(jday, rstep, lai, clumping, parameter, meteo, CosZs, var_o, var_n, p_soil, mid_res, mid_ET);
}
