/// @file snowpack.c
/// @brief This module will calculate the percentage of canopy and ground covered by snow
///        and output albedo of snow (used in energy balance) and density of snow in this step.
/// @details [snowpack_stage1] happens before any consumption of snow in this step, after the snow fall (supply)
/// @details [snowpack_stage2] happens after sublimation from ground and canopy (demand)
/// @details [snowpack stage3] happens after frozen and melt of snow pack (demand)
/// @author XZ Luo
/// @date May 25, 2015


# include "beps.h"


/// @brief Function of snowpack stage1.
/// @details [snowpack_stage1] happens before any consumption of snow in this step, after the snow fall (supply)
/// @details [Input] air temperature, precipitation,depth of snow from last step, density of snow from last step,
///                  mass of snow on canopy and ground (per ground area) from last step, length of step,
///                  leaf area index of overstorey and understorey excluding stem, albedo of snow from last step.
/// @details [Output] mass of snow on canopy and ground accumulation of snowfall,
///                   albedo of snow in this step,
///                   density of snow in this step.
/// @param  temp_air          air temperature
/// @param  precipitation     precipitation (m/s)
/// @param  mass_snow_o_last  weight of snow at overstorey from last step
/// @param  mass_snow_u_last  weight of snow at understorey from last step
/// @param  mass_snow_g_last  weight of snow on ground from last step
/// @param  mass_snow_o       mass of intercepted snow at overstory, input from last step, kg/m2
/// @param  mass_snow_u       mass of intercepted snow at understory, input from last step, kg/m2
/// @param  mass_snow_g       mass of intercepted snow on ground, input from last step, kg/m2
/// @param  lai_o             overstory lai
/// @param  lai_u             understory lai
/// @param  clumping          clumping index
/// @param  area_snow_o       area of snow at overstorey
/// @param  area_snow_u       area of snow at understorey
/// @param  percent_snow_o    percentage of snow cover at overstory, DECIDED by weight
/// @param  percent_snow_u    percentage of snow cover at understory, DECIDED by weight
/// @param  percent_snow_g    percentage of snow cover on ground, DECIDED by weight
/// @param  density_snow      density of snowpack on ground, input from last step, then changed in this module
/// @param  depth_snow        depth of snowpack, input from last step, changed here, then changed in stage2
/// @param  albedo_v_snow     visible albedo of snow, input from this step, changed in this module
/// @param  albedo_n_snow     near infrared albedo of snow, input from this step, changed in this module
/// @return void
void snowpack_stage1(double temp_air, double precipitation,double mass_snow_o_last, double mass_snow_u_last,
                     double mass_snow_g_last, double* mass_snow_o, double* mass_snow_u, double* mass_snow_g,
                     double lai_o, double lai_u, double clumping, double* area_snow_o, double* area_snow_u,
                     double* percent_snow_o, double* percent_snow_u, double* percent_snow_g,
                     double* density_snow, double* depth_snow, double* albedo_v_snow, double* albedo_n_snow)
{
    double massMax_snow_o, massMax_snow_u;     // maximum weight of snow on overstorey and understorey
    double massStep_snow_o, massStep_snow_u;   // weight (mass) of snow accumulated in this step

    double areaMax_snow_o, areaMax_snow_u;     // Maximum area of snow at overstorey and understorey

    double change_depth_snow;                  //change of snow depth on ground
    double density_water, density_new_snow;    // density of water, density of newly fallen snow
    double snowrate, snowrate_o, snowrate_u, snowrate_g; // snow rate, and snow rate for every part
    double albedo_v_Newsnow, albedo_n_Newsnow; // albedo of newly fallen snow in visible and near infrared band

    double length_step;

    length_step=kstep;

    density_new_snow=67.9+51.3*exp(temp_air/2.6);
    density_water=1025;
    albedo_v_Newsnow=0.94;
    albedo_n_Newsnow=0.8;

    massMax_snow_o=0.1*lai_o;
    massMax_snow_u=0.1*lai_u;
    areaMax_snow_o=lai_o*0.01;
    areaMax_snow_u=lai_u*0.01;

    // previous mass_snow 
    // mass_snow_o_last=*mass_snow_o;
    // mass_snow_u_last=*mass_snow_u;
    // mass_snow_g_last=*mass_snow_g;

    /*****************************/

    if (temp_air>0) // if temperature >0, no snow fall
        snowrate=0;
    else
        snowrate=precipitation*density_water/density_new_snow;

    snowrate_o = 0.0;

    if (temp_air<0) // if there is snow
    {
        //overstorey
        snowrate_o=snowrate;
        *mass_snow_o=mass_snow_o_last+snowrate_o*length_step*density_new_snow*(1-exp(-lai_o*clumping));

        *percent_snow_o = *mass_snow_o/massMax_snow_o;
        *percent_snow_o = max(0,*percent_snow_o);
        *percent_snow_o = min(1,*percent_snow_o);

        // change the weight based percentage to area based percentage

        *area_snow_o=*percent_snow_o*areaMax_snow_o;
        //*percentArea_snow_o=area_snow_o/(lai_o)/2;

        massStep_snow_o=*mass_snow_o-mass_snow_o_last;

        //understorey
        snowrate_u=snowrate_o-(massStep_snow_o)/density_new_snow/length_step;
        snowrate_u=max(0,snowrate_u);

        *mass_snow_u=mass_snow_u_last+snowrate_u*length_step*density_new_snow*(1-exp(-lai_u*clumping));

        *percent_snow_u = *mass_snow_u/massMax_snow_u;
        *percent_snow_u = max(0,*percent_snow_u);
        *percent_snow_u = min(1,*percent_snow_u);

        // change the weight based percentage to area based percentage

        *area_snow_u=*percent_snow_u*areaMax_snow_u;
        //*percentArea_snow_u=area_snow_u/(lai_u)/2;

        massStep_snow_u=*mass_snow_u-mass_snow_u_last;

        //ground, output mass of snow, percent of snow coverage, density of snow, depth of snow
        snowrate_g=snowrate_u-(massStep_snow_u)/density_new_snow/length_step;
        snowrate_g=max(0,snowrate_g);

        change_depth_snow=snowrate_g*length_step;
    }
    else
    {
        //overstorey
        *mass_snow_o = mass_snow_o_last;
        *percent_snow_o = *mass_snow_o/massMax_snow_o;
        *percent_snow_o = max(0,*percent_snow_o);
        *percent_snow_o = min(1,*percent_snow_o);

        *area_snow_o = *area_snow_o;

        // understorey
        *mass_snow_u = mass_snow_u_last;
        *percent_snow_u = *mass_snow_u/massMax_snow_u;
        *percent_snow_u = max(0,*percent_snow_u);
        *percent_snow_u = min(1,*percent_snow_u);
        *area_snow_u = *area_snow_u;

        //ground
        change_depth_snow=0;

    }
    change_depth_snow=max(0,change_depth_snow);
    *mass_snow_g=mass_snow_g_last+change_depth_snow*density_new_snow;
    *mass_snow_g=max(0,*mass_snow_g);

    if (change_depth_snow > 0)
        *density_snow=(*density_snow*(*depth_snow)+density_new_snow*change_depth_snow)/((*depth_snow)+change_depth_snow);
    else
        *density_snow=(*density_snow-250)*exp(-0.001*length_step/3600.0) + 250.0;

    if (*mass_snow_g > 0)
        *depth_snow=min(*mass_snow_g/(*density_snow), 10.0);
    else
        *depth_snow=0;

    *percent_snow_g=*mass_snow_g/(0.05*(*density_snow));
    *percent_snow_g=min(*percent_snow_g,1);

    // *deltaz=change_depth_snow;
    // albedo of snow in this step
    if (snowrate_o>0)
    {
        *albedo_v_snow=(*albedo_v_snow-0.70)*exp(-0.005*length_step/3600)+0.7;
        *albedo_n_snow=(*albedo_n_snow-0.42)*exp(-0.005*length_step/3600)+0.42;
    }
    else
    {
        *albedo_v_snow=albedo_v_Newsnow;
        *albedo_n_snow=albedo_n_Newsnow;
    }


}


/// @brief Function of snowpack stage2.
///        This module will calculate the snow remained on canopy surface after evaporation in this step
/// @details [snowpack_stage2] happens after sublimation from ground and canopy (demand)
/// @details [input] mass of snow on leaves after precipitation in this step,
///                  sublimation from leaves in this step
/// @details [output] mass of snow on leaves after the sublimation on leaves in this step
/// @param evapo_snow_o  evaporation of intercepted rain in this step, overstorey, kg/m2/s = mm/s
/// @param evapo_snow_u  evaporation of intercepted rain in this step, understorey, kg/m2/s = mm/s
/// @param mass_snow_o   supply of rain on leaves, overstorey, already added precipitation in this step
/// @param mass_snow_u   supply of rain on leaves, understorey, already added precipitation in this step
/// @return void
void snowpack_stage2(double evapo_snow_o, double evapo_snow_u,
                     double* mass_snow_o, double* mass_snow_u)

{
    double length_step; // length of step
    length_step=kstep;

    *mass_snow_o=*mass_snow_o-evapo_snow_o*length_step;
    *mass_snow_o=max(0,*mass_snow_o);

    *mass_snow_u=*mass_snow_u-evapo_snow_u*length_step;
    *mass_snow_u=max(0,*mass_snow_u);
}


/// @brief Function of snowpack stage3.
///        This module simulates the process of snow melting and water frozen in this step
/// @details [snowpack stage3] happens after frozen and melt of snow pack (demand)
/// @details [input] depth of snow on ground after stage 1, air temperature, ground surface temperature
/// @details [output] the amount of the melted snow, frozen snow
/// @param  temp_air        temperature of air in this step
/// @param  temp_snow       temperature of snow in this step
/// @param  temp_snow_last  temperature of snow in last step
/// @param  density_snow    density of snow output from stage1
/// @param  depth_snow      depth of snow on ground after stage1
/// @param  mass_snow_g     mass of snow on ground after stage1
/// @param  depth_water     depth of water after all precipitation and evaporation
/// @return void
void snowpack_stage3(double temp_air, double temp_snow, double temp_snow_last, double density_snow,
                     double* depth_snow, double* depth_water, double* mass_snow_g)
//double evapo_snow_g; // sublimation of snow on ground, input from evaporation_soil module, in mm/s
{
    double depth_snow_sup, mass_snow_sup; // depth and mass of snow after stage1, and minus the amount of sublimation. ( or the total supply of snow in this step)
    double mass_snow_melt, mass_water_frozen; // newly melted or frozen snow in this step, kg/m2;
    double cp_ice; // specific heat of ice 2228.261 J/kg/C
    double latent_fusion; // latent heat for fusion 3.34*1000000 J/kg
    double density_water;

    double melt_depth_snow, frozen_depth_snow; // change of snow weight valued in depth of snow
    double melt_depth_water, frozen_depth_water; // change of snow weight valued in depth of water

    double length_step; // length_step in model

    length_step=kstep;

    // it is assumed sublimation happens before the melting and freezing process
    depth_snow_sup=*depth_snow; // this snow depth has already considered sublimation
    mass_snow_sup=*mass_snow_g;

    //parameters
    cp_ice=2228.261;
    latent_fusion=3.34*1000000;
    density_water=1025;

    //simulate snow melt and freeze process
    mass_snow_melt=0;
    mass_water_frozen=0;

    // bool con_melt = temp_snow > 0 && temp_snow_last <= 0 && mass_snow_sup > 0;
    // bool con_froze = temp_snow <= 0 && temp_snow_last > 0 && *depth_water > 0;
    bool con_melt = temp_snow > 0 && mass_snow_sup > 0;
    bool con_froze = temp_snow <= 0 && *depth_water > 0;

    // case 1 depth of snow <0.02 m
    if (depth_snow_sup<=0.02)
    {
        if (temp_air>0 && depth_snow_sup >0)
        {
            mass_snow_melt=temp_air*0.0075*length_step/3600*0.3;
            mass_snow_melt=min(mass_snow_sup, mass_snow_melt); // the amount of melted snow could not be larger than supply
        }
        else
            mass_snow_melt=0;
    }

        //case 2 depth of snow > 0.02 < 0.05 m
    else if (depth_snow_sup > 0.02 && depth_snow_sup <=0.05)
    {
        if (con_melt) {
            mass_snow_melt=temp_snow*cp_ice*density_snow*depth_snow_sup/latent_fusion;
            mass_snow_melt=min(mass_snow_sup, mass_snow_melt); // the amount of melted snow could not be larger than supply
        } else
            mass_snow_melt=0;

        if (con_froze) {
            mass_water_frozen=(0-temp_snow)*cp_ice*density_snow*depth_snow_sup/latent_fusion;
            mass_water_frozen=min(mass_water_frozen,*depth_water*density_water);
        } else
            mass_water_frozen=0;
    }

        // case 3 depth of snow > 0.05 m
    else if(depth_snow_sup>0.05)
    {
        if (con_melt) {
            mass_snow_melt=temp_snow*cp_ice*density_snow*depth_snow_sup*0.02/latent_fusion;
            mass_snow_melt=min(mass_snow_sup, mass_snow_melt); // the amount of melted snow could not be larger than supply
        } else
            mass_snow_melt=0;

        if (con_froze) {
            mass_water_frozen=(0-temp_snow)*cp_ice*density_snow*depth_snow_sup*0.02/latent_fusion;
            mass_water_frozen=min(mass_water_frozen,*depth_water*density_water);
        } else
            mass_water_frozen=0;
    }

    // change in mass of snow on ground
    *mass_snow_g=*mass_snow_g-mass_snow_melt+mass_water_frozen;
    *mass_snow_g=max(0,*mass_snow_g);

    // change of depth in snow
    melt_depth_snow=mass_snow_melt/density_snow;
    frozen_depth_snow=mass_water_frozen/density_snow;
//    *depth_snow=depth_snow_sup-melt_depth_snow+frozen_depth_water;
    *depth_snow=depth_snow_sup-melt_depth_snow+frozen_depth_snow;
    *depth_snow=max(0,*depth_snow);

    // change of depth in water
    melt_depth_water=mass_snow_melt/density_water;
    frozen_depth_water=mass_water_frozen/density_water;
    *depth_water=*depth_water+melt_depth_water-frozen_depth_water;
    *depth_water=max(0,*depth_water);
}

