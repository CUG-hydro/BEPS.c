/// @file beps.h
/// @brief Header file for defining constants and global variables for BEPS program
///
/// CCRS (EMS/Applications Division)
///
/// @author Written by: J. Liu, Modified by:  G. Mo
/// @date June 2015

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "soil.h"

/// @brief Define Constants
#define NOERROR 0
#define ERROR 1
#define PI 3.1415926
#define zero 0.0000000001

// #define max(a,b) 	(a>b)?a:b         // used for UNIX
// #define min(a,b)	(a<b)?a:b	      // used for UNIX
#define max(a, b) ((a) > (b)) ? (a) : (b)  // LHE. the  original one can lead to disorder.
#define min(a, b) ((a) < (b)) ? (a) : (b)  // LHE

#define l_sta 105  // start line
#define l_end 105  // end line
#define p_sta 101  // start pix
#define p_end 101  // end pix

#define RTIMES 24    // 24
#define step 3600    // 3600 in sec
#define kstep 360    // 10 times per hour, 360 sec. per time
#define kloop 10     // 10 times per hour, 360 sec. per time
#define MAX_Loop 11  // 10 times per hour, 360 sec. per time
#define layer 5
#define depth_f 6
#define CO2_air 380  // atmospheric CO2 concentration
#define rho_a 1.292  // density of air at 0C

/// @brief Declare structures
struct climatedata {
    double Srad;
    double LR;
    double temp;
    double rh;
    double rain;
    double wind;
    double dr_o;
    double df_o;
    double dr_u;
    double df_u;
    //	float st_c;
};

struct results {
    double gpp_o_sunlit;
    double gpp_u_sunlit;
    double gpp_o_shaded;
    double gpp_u_shaded;
    double plant_resp;
    double npp_o;
    double npp_u;
    double GPP;
    double NPP;
    double NEP;
    double soil_resp;
    double Net_Rad;
    double SH;
    double LH;
    double Trans;
    double Evap;
};

struct OutputET {
    double Trans_o;
    double Trans_u;
    double Eil_o;
    double Eil_u;
    double EiS_o;
    double EiS_u;
    double Evap_soil;
    double Evap_SW;
    double Evap_SS;
    double Qhc_o;
    double Qhc_u;
    double Qhg;
    double Trans;
    double Evap;
    double SH;
    double LH;
};

struct cpools {
    double Ccd[3];
    double Cssd[3];
    double Csmd[3];
    double Cfsd[3];
    double Cfmd[3];
    double Csm[3];
    double Cm[3];
    double Cs[3];
    double Cp[3];
};

/** Leaf structure ----------------------------------------------------------**/
typedef struct {
    double o_sunlit;
    double o_shaded;
    double u_sunlit;
    double u_shaded;
} Leaf;

// Leaf Leaf_default = {0, 0, 0, 0};
// #define LEAF(...) ((Leaf){.x = 0, .y = 0, ##__VA_ARGS__})

void init_leaf_struct(Leaf* x, Leaf replacement);
void init_leaf_dbl(Leaf* x, double replacement);
void init_leaf_dbl2(Leaf* x, double overstory, double understory);

Leaf leaf_mul(Leaf x, Leaf y);
double leaf_sum(Leaf x);

/// @brief Declare functions
// void readconf();
// void mid_prg();
// void readinput1();
// void readlai_d();
// void readlonlat();

void update_ET(struct OutputET* x, struct results* mid_res, double Ta);

void inter_prg_c(int jday, int rstep, double lai, double clumping, double parameter[], struct climatedata* meteo,
                                                                                        double CosZs, double var_o[], double var_n[], struct Soil* soilp, struct results* mid_res, struct OutputET* mid_ET);

void inter_prg(int jday, int rstep, double lai, double clumping, double parameter[], struct climatedata* meteo,
               double CosZs, double var_o[], double var_n[], struct Soil* soilp, struct results* mid_res);

void s_coszs(short jday, short j, float lat, float lon, double* CosZs);

void aerodynamic_conductance(double canopy_height_o, double canopy_height_u, double zz, double clumping,
                             double temp_air, double wind_sp, double SH_o_p, double lai_o, double lai_u,
                             double* rm, double* ra_u, double* ra_g, double* G_o_a, double* G_o_b, double* G_u_a, double* G_u_b);
void plantresp(int LC, struct results* mid_res, double lai_yr, double lai, double temp_air, double temp_soil, double CosZs);

void Vcmax_Jmax(double lai_o, double clumping, double Vcmax0,
                double slope_Vcmax_N, double leaf_N, double CosZs,
                double* Vcmax_sunlit, double* Vcmax_shaded, double* Jmax_sunlit, double* Jmax_shaded);

void soilresp(double* Ccd, double* Cssd, double* Csmd, double* Cfsd, double* Cfmd,
              double* Csm, double* Cm, double* Cs, double* Cp, float npp_yr, double* coef,
              int soiltype, struct Soil* soilp, struct results* mid_res);

void readparam(short lc, double parameter1[]);

void lai2(double clumping, double CosZs, double stem_o, double stem_u, double lai_o, double lai_u,
          Leaf* LAI, Leaf* PAI);

void readcoef(short lc, int stxt, double coef[]);

void readhydr_param();

void photosynthesis(double temp_leaf_p, double rad_leaf, double e_air, double g_lb_w, double vc_opt,
                    double f_soilwater, double b_h2o, double m_h2o, double cii, double temp_leaf_c, double LH_leaf,
                    double* Gs_w, double* aphoto, double* ci);
void soil_water_factor();

/** using Leaf struct */
void netRadiation(double shortRad_global, double CosZs, double temp_o, double temp_u, double temp_g,
                  double lai_o, double lai_u, double lai_os, double lai_us,  // LAI of overstorey and understorey, with and without stem
                  Leaf lai,
                  double clumping, double temp_air, double rh,
                  double albedo_snow_v, double albedo_snow_n, double percentArea_snow_o, double percentArea_snow_u, double percent_snow_g,
                  double albedo_v_o, double albedo_n_o, double albedo_v_u, double albedo_n_u, double albedo_v_g, double albedo_n_g,
                  double* netRad_o, double* netRad_u, double* netRad_g,
                  Leaf* netRadLeaf,
                  Leaf* netShortRadLeaf);

void Leaf_Temperatures(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca,
                       Leaf Gw, Leaf Gww, Leaf Gh,
                       double Xcs_o, double Xcl_o, double Xcs_u, double Xcl_u,
                       Leaf radiation, Leaf* Tc);
double Leaf_Temperature(double Tair, double slope, double psychrometer, double VPD_air, double Cp_ca,
                        double Gw, double Gww, double Gh, double Xc_sl, double radiation, bool constrain);

void sensible_heat(
    Leaf tempL, double temp_g, double temp_air, double rh_air,
    Leaf Gheat, double Gheat_g, Leaf lai,
    double* SH_o, double* SH_u, double* SH_g);

void transpiration(Leaf tempL, double temp_air, double rh_air,
                   Leaf Gtrans, Leaf lai,
                   double* trans_o, double* trans_u);

void evaporation_canopy(
    Leaf tempL, double temp_air, double rh_air,
    Leaf Gwater, Leaf lai,
    double percent_water_o, double percent_water_u, double percent_snow_o, double percent_snow_u,
    double* evapo_water_o, double* evapo_water_u, double* evapo_snow_o, double* evapo_snow_u);

/** end of Leaf struct usage */

void evaporation_soil(double temp_air, double temp_g, double rh_air, double netRad_g, double Gheat_g,
                      double* percent_snow_g, double* depth_water, double* depth_snow, double* mass_water_g, double* mass_snow_g,
                      double density_snow, double swc_g, double porosity_g,
                      double* evapo_soil, double* evapo_water_g, double* evapo_snow_g);

void rainfall_stage1(double temp_air, double precipitation, double mass_water_o_last, double mass_water_u_last,
                     double lai_o, double lai_u, double clumping,
                     double* mass_water_o, double* mass_water_u,
                     double* percent_water_o, double* percent_water_u, double* precipitation_g);

void rainfall_stage2(double evapo_water_o, double evapo_water_u,
                     double* mass_water_o, double* mass_water_u);
void rainfall_stage3();

void meteo_pack(double temp, double rh, double* meteo_pack_output);

void surface_temperature(double temp_air, double rh_air, double depth_snow, double depth_water,
                         double capacity_heat_soil1, double capacity_heat_soil0, double Gheat_g,
                         double depth_soil1, double density_snow, double tempL_u, double netRad_g,
                         double evapo_soil, double evapo_water_g, double evapo_snow_g, double lambda_soil1,
                         double percent_snow_g, double heat_flux_soil1, double temp_ground_last,
                         double temp_soil1_last, double temp_any0_last, double temp_snow_last,
                         double temp_soil0_last, double temp_snow1_last, double temp_snow2_last,
                         double* temp_ground, double* temp_any0, double* temp_snow,
                         double* temp_soil0, double* temp_snow1, double* temp_snow2, double* heat_flux);

void snowpack_stage1(double temp_air, double precipitation, double mass_snow_o_last, double mass_snow_u_last,
                     double mass_snow_g_last, double* mass_snow_o, double* mass_snow_u, double* mass_snow_g,
                     double lai_o, double lai_u, double clumping, double* area_snow_o, double* area_snow_u,
                     double* percent_snow_o, double* percent_snow_u, double* percent_snow_g,
                     double* density_snow, double* depth_snow, double* albedo_v_snow, double* albedo_n_snow);

void snowpack_stage2(double evapo_snow_o, double evapo_snow_u,
                     double* mass_snow_o, double* mass_snow_u);

void snowpack_stage3(double temp_air, double temp_snow, double temp_snow_last, double density_snow,
                     double* depth_snow, double* depth_water, double* mass_snow_g);

// /// @brief Declare global variables
// short lc_no;
// int yr,bgn_day,end_day;
// int npixels,nlines;

// char lc_fn[255];       /* Land cover file */
// char lai_fn[255];      /* Leaf area index file  */
// char lai_fp[255];      /* Leaf area index file prefix */
// char stxt_fn[255];     /* soil texture file */
// char ci_fn[255];       /* clumping index file */
// char st_fn[255];       /* initial values of soil temp */
// char sw_fn[255];       /* initial values of soil water */
// char sdp_fn[255];      /* initial values of snow depth*/

// char r_fn[255];        /* meteor. data files */
// char t_fn[255];
// char h_fn[255];
// char p_fn[255];
// char wd_fn[255];

// char lon_fn[255];
// char lat_fn[255];

// char fp4outp1[255];     /* output file1 prefix */
// char fp4outp2[255];     /* output file2 prefix */
// char fp4outp3[255];     /* output file3 prefix */
