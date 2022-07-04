import Parameters: @with_kw, @with_kw_noshow


zero10() = tuple(zeros(10)...)

@with_kw mutable struct Soil
    flag::Cint = Cint(0)
    n_layer::Cint = Cint(5)
    step_period::Cint  = Cint(1)
    Zp::Cdouble = Cdouble(0)
    Zsp::Cdouble = Cdouble(0)
    r_rain_g::Cdouble = Cdouble(0)
    soil_r::Cdouble = Cdouble(0)
    r_drainage::Cdouble = Cdouble(0)
    r_root_decay::Cdouble = Cdouble(0)
    psi_min::Cdouble = Cdouble(0)
    alpha::Cdouble = Cdouble(0)
    f_soilwater::Cdouble = Cdouble(0)
    d_soil::NTuple{10, Cdouble} = zero10()
    f_root::NTuple{10, Cdouble} = zero10()
    dt::NTuple{10, Cdouble} = zero10()
    thermal_cond::NTuple{10, Cdouble} = zero10()
    theta_vfc::NTuple{10, Cdouble} = zero10()
    theta_vwp::NTuple{10, Cdouble} = zero10()
    fei::NTuple{10, Cdouble} = zero10()
    Ksat::NTuple{10, Cdouble} = zero10()
    psi_sat::NTuple{10, Cdouble} = zero10()
    b::NTuple{10, Cdouble} = zero10()
    density_soil::NTuple{10, Cdouble} = zero10()
    f_org::NTuple{10, Cdouble} = zero10()
    ice_ratio::NTuple{10, Cdouble} = zero10()
    thetam::NTuple{10, Cdouble} = zero10()
    thetam_prev::NTuple{10, Cdouble} = zero10()
    temp_soil_p::NTuple{10, Cdouble} = zero10()
    temp_soil_c::NTuple{10, Cdouble} = zero10()
    f_ice::NTuple{10, Cdouble} = zero10()
    psim::NTuple{10, Cdouble} = zero10()
    thetab::NTuple{10, Cdouble} = zero10()
    psib::NTuple{10, Cdouble} = zero10()
    r_waterflow::NTuple{10, Cdouble} = zero10()
    km::NTuple{10, Cdouble} = zero10()
    Kb::NTuple{10, Cdouble} = zero10()
    KK::NTuple{10, Cdouble} = zero10()
    Cs::NTuple{10, Cdouble} = zero10()
    lambda::NTuple{10, Cdouble} = zero10()
    Ett::NTuple{10, Cdouble} = zero10()
    G::NTuple{10, Cdouble} = zero10()
end

function SoilRootFraction(soil)
    ccall((:SoilRootFraction, libbeps), Cvoid, (Ptr{Soil},), soil)
end

function Init_Soil_Parameters(landcover::Int, stxt::Int, r_root_decay::Real, p::Soil)
    ccall((:Init_Soil_Parameters, libbeps), Cvoid, (Cint, Cint, Cdouble, Ptr{Soil}), 
    Cint(landcover), Cint(stxt), Cfloat(r_root_decay), Ref(p))
end

function Init_Soil_Status(p, Tsoil, Tair, Ms, snowdepth)
    ccall((:Init_Soil_Status, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble), p, Tsoil, Tair, Ms, snowdepth)
end

function soil_water_factor_v2(p)
    ccall((:soil_water_factor_v2, libbeps), Cvoid, (Ptr{Soil},), p)
end

function Soil_Water_Uptake(p, Trans_o, Trans_u, Evap_soil)
    ccall((:Soil_Water_Uptake, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble), p, Trans_o, Trans_u, Evap_soil)
end

function UpdateSoilLambda(soil)
    ccall((:UpdateSoilLambda, libbeps), Cvoid, (Ptr{Soil},), soil)
end

function init_soil_parameter(T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, soil)
    ccall((:init_soil_parameter, libbeps), Cvoid, (Cuchar, Cuchar, Cuchar, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Soil}), T_USDA, S_USDA, Ref_Depth, T_Density, S_Density, T_OC, S_OC, soil)
end

function Update_Cs(p)
    ccall((:Update_Cs, libbeps), Cvoid, (Ptr{Soil},), p)
end

function Update_ice_ratio(p)
    ccall((:Update_ice_ratio, libbeps), Cvoid, (Ptr{Soil},), p)
end

function UpdateSoilThermalConductivity(p)
    ccall((:UpdateSoilThermalConductivity, libbeps), Cvoid, (Ptr{Soil},), p)
end

function UpdateHeatFlux(p, Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
    ccall((:UpdateHeatFlux, libbeps), Cvoid, (Ptr{Soil}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), p, Xg_snow, lambda_snow, Tsn0, Tair_annual_mean, peroid_in_seconds)
end

function UpdateSoilMoisture(p, peroid_in_seconds)
    ccall((:UpdateSoilMoisture, libbeps), Cvoid, (Ptr{Soil}, Cdouble), p, peroid_in_seconds)
end
