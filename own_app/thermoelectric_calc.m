%% Load essential parameters
% Implement struct data structure in the future!

% warning('off','all');           % Turn off all warnings
run("param_thermoelectric_cooling.m");

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global Area_cross_sect_cold Dh_cold
global R_e_hc R_k_hc alpha_seeback num_semi_cond
global fin_width_cold fin_length_cold fin_thickness_cold sink_height_cold num_fins_cold k_fin_cold per_fin_area_cold base_area_cold fin_area_total_cold
global fin_width_hot fin_length_hot fin_thickness_hot sink_height_hot num_fins_hot k_fin_hot per_fin_area_hot base_area_hot fin_area_total_hot 

%% Define simulation parameters (CHANGME)

% General parameters
J_e = 0;              % Optimal current (CHANGE TO FUNCTION)
J_iters = 20;
J_max = 4.0;

% Initial conditions - Cold Side (Air restricted to channel)
inlet_temp_cold = 308.15;   % K
air_speed_cold = 2.2;      % m/s
m_dot_air_cold = Area_cross_sect_cold * rho_air * air_speed_cold;

% Initial conditions - Hot Side (Air not restricted to channel)
inlet_temp_hot = 308.15;   % K
CFM_fan_hot = 48;              % CubicFt/min
volumetric_flow_rate_hot = CFM_fan_hot * ((0.3048^3) / 60);   % m^3/s - conversion factor
m_dot_air_hot = volumetric_flow_rate_hot / rho_air;
fan_area = pi * 0.04^2;
air_speed_hot = volumetric_flow_rate_hot / fan_area ;         % m/s

% Compute convective coefficient & fin efficiencies
[R_ku_cold, h_cold] = compute_convective_coefficient_cold_NTU(air_speed_cold, fin_area_total_cold, fin_width_cold, Dh_cold, m_dot_air_cold);
[R_ku_hot, h_hot] = compute_convective_coefficient_hot_without_NTU(air_speed_hot, fin_area_total_hot, fin_width_hot);
overall_fin_eff_hot = compute_fin_efficiency(h_hot, k_fin_hot, fin_thickness_hot, fin_length_hot, num_fins_hot, per_fin_area_hot, fin_area_total_hot);        
overall_fin_eff_cold = compute_fin_efficiency(h_cold, k_fin_cold, fin_thickness_cold, fin_length_cold, num_fins_cold, per_fin_area_cold, fin_area_total_cold);        


%% Print initialization message
fprintf('<strong>***Initialization***\n</strong>');
fprintf('Inlet Air Temperature - Cold Side (T_in_cold): %.3f K \n', inlet_temp_cold);
fprintf('Inlet Air Speed - Cold Side (U_cold): %.1f m/s \n', air_speed_cold);
fprintf('Inlet Air Temperature - Hot Side (T_in_hot): %.3f K \n', inlet_temp_hot);
fprintf('Inlet Air Speed - Hot Side (U_hot): %.1f m/s \n', air_speed_hot);
fprintf('Convective Coefficient Resistance (R_ku_c) - Cold Side: %.3f K/W\n', R_ku_cold);
fprintf('Convective Coefficient Resistance (R_ku_h) - Hot Side: %.3f K/W\n', R_ku_hot);
fprintf('Conductive Coefficient Resistance (R_k_hc): %.3f K/W\n\n', R_k_hc);

%% Main Calculation Body

cooling_power_arr = zeros(J_iters, 1);
power_required_arr = zeros(J_iters, 1);
delta_J_arr = linspace(0, J_max, J_iters);
J_optimal = 0;
max_cooling_power = 0;
T_h_optimal = 0;
T_c_optimal = 0;
power_required_optimal = 0;
COP_optimal = 0;
outlet_temp_cold_optimal = 0;
outlet_temp_hot_optimal = 0;

for i = 1:length(delta_J_arr)
    
    J_e = delta_J_arr(i);
    
    % x = T_h, y = T_c, z = Q_c
    syms x y z
    eqn1 = ((x - y) / R_k_hc) + (overall_fin_eff_hot * (x-inlet_temp_hot) / R_ku_hot) == (num_semi_cond * alpha_seeback * J_e * x) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
    eqn2 = (-(x - y) / R_k_hc) + z == (-num_semi_cond * alpha_seeback * J_e * y) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
    eqn3 = z == (overall_fin_eff_cold * (y - inlet_temp_cold) ) / R_ku_cold;

    sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
    T_h_peltier = double(sol.x);
    T_c_peltier = double(sol.y);
    Q_c_peltier = double(sol.z);
    
    Q_h_peltier = (inlet_temp_hot - T_h_peltier) / R_ku_hot;
    power_conduction_peltier = (T_h_peltier - T_c_peltier) / R_k_hc;
    outlet_temp_cold = inlet_temp_cold + Q_c_peltier/(m_dot_air_cold * Cp_air);
    outlet_temp_hot = inlet_temp_hot - Q_h_peltier/(m_dot_air_hot * Cp_air);
    power_required = num_semi_cond * ((R_e_hc * J_e^2) + (alpha_seeback * J_e * (T_h_peltier - T_c_peltier)) );
    coefficient_performance = -100 * Q_c_peltier / power_required;
    
    cooling_power_arr(i) = Q_c_peltier;
    power_required_arr(i) = power_required;
    
    % Find optimal current which gives max cooling
    if -Q_c_peltier > -max_cooling_power
        max_cooling_power = Q_c_peltier;
        max_heating_power = Q_h_peltier;
        J_optimal = J_e;
        T_h_optimal = T_h_peltier;
        T_c_optimal = T_c_peltier;
        power_required_optimal = power_required;
        COP_optimal = coefficient_performance;
        outlet_temp_cold_optimal = outlet_temp_cold;
        outlet_temp_hot_optimal = outlet_temp_hot;
    end
    
    % Print results..
    fprintf('<strong>===Iteration %d===\n</strong>', i);
    fprintf('Input Current (J_e): %.1f A \n', J_e);
    fprintf('Power Required (P_e): %.1f W\n', power_required);
    fprintf('Coefficient of Performance (COP): %.1f %% \n', coefficient_performance);
    
    % Cold side
    fprintf('Cold side Temperature (T_c): %.1f K \n', T_c_peltier);
    fprintf('Cooling Power - Cold Side (Q_c_peltier): %.2f W\n', Q_c_peltier);
    fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.1f K\n', outlet_temp_cold);
    
    % Hot side
    fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_peltier);
    fprintf('Heating Power - Hot Side (Q_h_peltier): %.2f W\n', Q_h_peltier);
    fprintf('Outlet Air Temperature - Hot Side (T_out_hot): %.1f K\n\n', outlet_temp_hot);

end

% Display optimal results
% outlet_temp_cold_optimal = inlet_temp_cold + max_cooling_power/(m_dot_air_cold * Cp_air);

fprintf('<strong>===Optimal Results===\n</strong>');

fprintf('Input Current (J_e): %.1f A \n', J_optimal);
fprintf('Power Required (P_e): %.1f W\n', power_required_optimal);
fprintf('Coefficient of Performance (COP): %.1f %% \n', COP_optimal);

% Cold side
fprintf('Cold side Temperature (T_c): %.1f K \n', T_c_optimal);
fprintf('Cooling Power - Cold Side (Q_c_peltier): %.2f W\n', max_cooling_power);
fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.1f K\n', outlet_temp_cold_optimal);

% Hot side
fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_optimal);
fprintf('Heating Power - Hot Side (Q_h_peltier): %.2f W\n', max_heating_power);
fprintf('Outlet Air Temperature - Hot Side (T_out_hot): %.1f K\n\n', outlet_temp_hot_optimal);


%% Plot final graphs

% Plot graph of Cooling power against Current
figure(1)
plot(delta_J_arr, cooling_power_arr);
title("Cooling Power against Input Current");
xlabel("Current [A]");
ylabel("Cooling Power [W]");
grid on;

% Plot abs cooling power and power consumption against current
hold on;
figure(2)
plot(delta_J_arr, -cooling_power_arr, delta_J_arr, power_required_arr);
title("Absolute Cooling Power and Power Consumption against Input Current");
xlabel("Current [A]");
ylabel("Power [W]");
legend("Cooling Power", "Power Consumed", "Location", "NorthEast");
grid on;

%% Main Functions Used


% Assuming flow over plate (Likely laminar Re < 5 * 10^5)
function [R_ku, h] = compute_convective_coefficient_cold_NTU(air_speed, Area_fin_total, fin_width, Dh, mass_dot)
    
    global kin_visc_air k_air Pr_air Cp_air;

    Re = (air_speed * fin_width)/kin_visc_air;
    Nu = 0.664 * Re^(0.5) * Pr_air^(1/3);
    NTU = (Area_fin_total * Nu * k_air) / (fin_width * mass_dot * Cp_air);
    heat_transfer_eff = 1 - exp(-NTU);                      % R changes with flow in channel...
    R_ku = 1 / (mass_dot * Cp_air * heat_transfer_eff);     % Average convective resistance
    h = 1 / (R_ku * Area_fin_total);                        % R = 1/hA
    
%     R_ku = fin_width/(Area_fin_total * Nu * k_air); 
%     h = Nu * k_air / fin_width;

end

% Assuming flow over plate (Likely laminar Re < 5 * 10^5)
function [R_ku, h] = compute_convective_coefficient_hot_without_NTU(air_speed, Area_fin_total, fin_width)
    
    global kin_visc_air k_air Pr_air;

    Re = (air_speed * fin_width)/kin_visc_air;
    Nu = 0.664 * Re^(0.5) * Pr_air^(1/3);    
    R_ku = fin_width/(Area_fin_total * Nu * k_air); 
    h = Nu * k_air / fin_width;

end

% Returns fin efficiency overall (including base and fins)
% Convection Coeff = h (Without considering area)
% conduction_coeff = k
function fin_eff_overall = compute_fin_efficiency(convection_coeff, conduction_coeff, fin_thickness, fin_length, num_fins, area_fin, area_total)
    
    m = sqrt(2 * convection_coeff / (conduction_coeff * fin_thickness) );
    L_c = fin_length + fin_thickness/2;
    fin_eff = (tanh(m*L_c) ) / (m * L_c);
    
    fin_eff_overall = 1 - ( (num_fins * area_fin / area_total) * (1 - fin_eff) );

end



