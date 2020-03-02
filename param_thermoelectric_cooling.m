%% Properties for calculation

clear global;
clc;clear;

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global height width Area_cross_sect perimeter Dh
global R_e_hc R_k_hc alpha_seeback num_semi_cond

%% Properties

% Table C.22 (300K) - Air properties
kin_visc_air = 1.566 * 10^-5;        % Calculate Reynolds (m^2/s)
Cp_air = 1005;                       % heat capacity (J/kgK)
k_air = 0.0267;                      % To calculate Conv Coeff (W/mK)
alpha_air = 2.257 * 10^-5;           % To calculate Le (m^2/s)
Pr_air = 0.69;                       % Prandt Num ()
rho_air = 1.177;                     % Density (Kg/m^3)


% Table C.9 - Bismuth Telluride Peltier element properties
alpha_s_pos = 2.3 * 10^-4;           % V/degC
alpha_s_neg = -2.1 * 10^-4;           % V/degC
rho_e_pos = 10^-5;                   % Ohm-m
rho_e_neg = rho_e_pos;               % Ohm-m
k_bismuth_pos = 1.7;                 % W/mK
k_bismuth_neg = 1.45;                % W/mK
width_semi_cond = 0.0015;            % m (Assume square base)
height_semi_cond = 0.0035;           % m
num_semi_cond = 400;

alpha_seeback = alpha_s_pos - alpha_s_neg;
R_e_hc = (height_semi_cond/(width_semi_cond^2)) * (rho_e_pos + rho_e_neg);
R_k_hc = 1 / ( num_semi_cond * (width_semi_cond^2/height_semi_cond) * (k_bismuth_pos + k_bismuth_neg) );

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