%% Load essential parameters

% warning('off','all');           % Turn off all warnings
run("param_thermoelectric_cooling.m");

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global height width Area_cross_sect perimeter Dh
global R_e_hc R_k_hc alpha_seeback num_semi_cond

%% Define simulation parameters (CHANGME)

% General parameters
J_e = 0.8;              % Optimal current (CHANGE TO FUNCTION)

% Initial conditions - Cold Side 
inlet_temp_cold = 293.15;   % K
air_speed_cold = 1;      % m/s
m_dot_air_cold = Area_cross_sect * rho_air * air_speed_cold;

% Initial conditions - Hot Side 
inlet_temp_hot = 298.15;   % K
air_speed_hot = 1;      % m/s  

% Fin conditions - Cold Side
fin_area_total_cold = 0.05;      % Given by prof's example [m]       
fin_width_cold = 0.09;           % length parallel to flow [m]
overall_fin_eff_cold = 1;        % (CHANGE TO FUNCTION)

% Fin conditions - Hot Side (ASSUMING 2* BIGGER ON ALL SIDES)
fin_area_total_hot = fin_area_total_cold*4;      % Given by prof's example [m]       
fin_width_hot = fin_width_cold*2;                % length parallel to flow [m]
overall_fin_eff_hot = 1;                         % (CHANGE TO FUNCTION)


%% Compute convective coefficient (resistance)

R_ku_cold = compute_convective_coefficient(air_speed_cold, fin_area_total_cold, fin_width_cold);
R_ku_hot = compute_convective_coefficient(air_speed_hot, fin_area_total_hot, fin_width_hot);
fprintf('<strong>***Initialization***\n</strong>');
fprintf('Inlet Air Temperature (T_in): %.3f K \n', inlet_temp_cold);
fprintf('Inlet Air Speed (U): %.1f m/s \n', air_speed_cold);
fprintf('Convective Coefficient Resistance (R_ku_c) - Cold Side: %.3f K/W\n', R_ku_cold);
fprintf('Convective Coefficient Resistance (R_ku_h) - Hot Side: %.3f K/W\n', R_ku_hot);
fprintf('Conductive Coefficient Resistance (R_k_hc): %.3f K/W\n\n', R_k_hc);

%% Main Calculation Body

% x = T_h, y = T_c, z = Q_c, 
syms x y z
eqn1 = ((x - y) / R_k_hc) + ((x-inlet_temp_hot)/1.515) == (num_semi_cond * alpha_seeback * J_e * x) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
eqn2 = (-(x - y) / R_k_hc) + z == (-num_semi_cond * alpha_seeback * J_e * y) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
eqn3 = z == (y - inlet_temp_cold) / 0.4183;

sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
T_h_peltier = double(sol.x);
T_c_peltier = double(sol.y);
Q_c_peltier = double(sol.z);

outlet_temp_cold = inlet_temp_cold + Q_c_peltier/(m_dot_air_cold * Cp_air);
coefficient_performance = -100 * Q_c_peltier / ( 400 * ((R_e_hc * J_e^2) + (alpha_seeback * J_e * (T_h_peltier - T_c_peltier)) ) );


%% Print results..
fprintf('<strong>===RESULTS===\n</strong>');
fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_peltier);
fprintf('Cold side Temperature (T_h): %.1f K \n', T_c_peltier);
fprintf('Cooling Power (Q_c_peltier): %.2f W\n', Q_c_peltier);
fprintf('Outlet Air Temperature (T_out): %.1f K\n', outlet_temp_cold);
fprintf('Coefficient of Performance (COP): %.1f %% \n\n', coefficient_performance);





%% Main Functions Used





% Assuming flow over plate (Likely laminar Re < 5 * 10^5)
function R_ku = compute_convective_coefficient(air_speed, Area_fin_total, fin_width)
    
    global kin_visc_air k_air Pr_air;

    Re = (air_speed * fin_width)/kin_visc_air;
    Nu = 0.664 * Re^(0.5) * Pr_air^(1/3);
    R_ku = fin_width/(Area_fin_total * Nu * k_air); 
%     Nu = 0.023 * Re^(4/5) * Pr_air^(0.3);       % n = 0.3
%     R_ku = Dh/(Area_water_contact * Nu * k_air);     

end



