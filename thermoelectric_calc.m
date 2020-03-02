%% Load essential parameters

warning('off','all');           % Turn off all warnings
run("param_thermoelectric_cooling.m");

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global height width Area_cross_sect perimeter Dh
global R_e_hc R_k_hc

%% Define simulation parameters (move some of these to param_evap_cool.m)

% Initial conditions - Cold Side (CHANGME)
inlet_temp_cold = 308;   % K
air_speed_cold = 1;      % m/s

% Initial conditions - Hot Side (CHANGME)
inlet_temp_hot = 308;   % K
air_speed_hot = 1;      % m/s  

% Fin conditions - Cold Side
fin_area_total_cold = 0.05;      % Given by prof's example [m]       
fin_width_cold = 0.09;           % length parallel to flow [m]

% Fin conditions - Hot Side (ASSUMING 2* BIGGER ON ALL SIDES)
fin_area_total_hot = fin_area_total_cold*4;      % Given by prof's example [m]       
fin_width_hot = fin_width_cold*2;           % length parallel to flow [m]



%% Compute convective coefficient (resistance)

R_ku_cold = compute_convective_coefficient(air_speed_cold, fin_area_total_cold, fin_width_cold);
R_ku_hot = compute_convective_coefficient(air_speed_hot, fin_area_total_hot, fin_width_hot);
fprintf('<strong>***Initialization***\n</strong>');
fprintf('Inlet Air Temperature (T_in): %.3f K \n', inlet_temp_cold);
fprintf('Inlet Air Speed (U): %.1f m/s \n', air_speed_cold);
fprintf('Convective Coefficient Resistance (R_ku) - Cold Side: %.3f K/W\n', R_ku_cold);
fprintf('Convective Coefficient Resistance (R_ku) - Hot Side: %.3f K/W\n', R_ku_hot);
fprintf('Conductive Coefficient Resistance (R_k_hc) - Hot Side: %.3f K/W\n\n', R_k_hc);

%% Main Calculation Body



%% Main Functions Used





% Assuming flow over plate (Likely laminar Re < 5 * 10^5)
function R_ku = compute_convective_coefficient(air_speed, Area_fin_total, fin_width)
    
    global kin_visc_air k_air Pr_air;

    Re = (air_speed * fin_width)/kin_visc_air;
    Nu = 0.664 * Re^(0.5) * Pr_air^(1/3);       % n = 0.3
    R_ku = fin_width/(Area_fin_total * Nu * k_air); 
%     Nu = 0.023 * Re^(4/5) * Pr_air^(0.3);       % n = 0.3
%     R_ku = Dh/(Area_water_contact * Nu * k_air);     

end



