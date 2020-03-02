%% Properties for calculation

clear global;
clc;clear;

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air T_lg delta_h_lg
global rho_water_in_air_ref height width Area_cross_sect perimeter Dh
global mol_weight_water R_g

%% Properties

% Table C.22 (300K)
kin_visc_air = 1.566 * 10^-5;            % Calculate Reynolds (m^2/s)
Cp_air = 1005;                       % heat capacity (J/kgK)
k_air = 0.0267;                      % To calculate Conv Coeff (W/mK)
alpha_air = 2.257 * 10^-5;           % To calculate Le (m^2/s)
Pr_air = 0.69;                           % Prandt Num ()
rho_air = 1.177;                     % Density (Kg/m^3)




%% Dimensions of flow channel (between sheets)
% Box dimension of 0.2m * 0.2m * 0.2m (evaCooler)

width = 0.12;
height = 0.12;
Area_cross_sect = height * width;
perimeter = (2 * height) + (2 * width);

Dh = 4*Area_cross_sect/perimeter; 

%% Peltier thermoelectric cooling

% mass_flow_air = 1 * (0.12^2) * rho_air;
% outlet_air_temp = 308 - (85 / (mass_flow_air * Cp_air) )
% deltaT = 308 - outlet_air_temp