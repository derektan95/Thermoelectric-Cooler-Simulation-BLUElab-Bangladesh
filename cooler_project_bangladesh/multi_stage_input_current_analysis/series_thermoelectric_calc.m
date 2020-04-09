%% Load essential parameters
% Implement struct data structure in the future!

% warning('off','all');           % Turn off all warnings
run("series_param_thermoelectric_cooling.m");

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global Area_cross_sect_cold_per_channel Dh_cold_per_channel num_channels area_per_channel
global R_e_hc R_k_hc alpha_seeback num_semi_cond
global fin_width_cold fin_length_cold fin_thickness_cold sink_height_cold num_fins_cold k_fin_cold per_fin_area_cold base_area_cold fin_area_total_cold
global fin_width_hot fin_length_hot fin_thickness_hot sink_height_hot num_fins_hot k_fin_hot per_fin_area_hot base_area_hot fin_area_total_hot 

%% Define simulation parameters (CHANGME)

% Number of stages (for peltier modules in series)
num_stages = 2;

% General parameters
J_e = 0;              % Optimal current (CHANGE TO FUNCTION)
J_iters = 100;
J_max = 10.0;

% Initial conditions - Cold Side (Air restricted to channel)
inlet_temp_cold = 308.15;   % K
CFM_nominal_cold = 33;                           % Nominal from specsheet (Small Fan = 5.8579, Large Fan = 59)
input_voltage_adjust_factor = 1;                                % Divide CFM by voltage divident
CFM_fan_cold = CFM_nominal_cold / input_voltage_adjust_factor;       % CubicFt/min (CFM_max = 5.8579)
volumetric_flow_rate_cold = CFM_fan_cold * ((0.3048^3) / 60);   % m^3/s - conversion factor
m_dot_air_cold = volumetric_flow_rate_cold * rho_air;
% fan_area_cold = pi * 0.02^2;                                   % CHANGEME
total_cross_section_fin_area = Area_cross_sect_cold_per_channel * num_channels;
air_speed_cold = volumetric_flow_rate_cold / total_cross_section_fin_area ;         % m/s

% air_speed_cold = 2.2;      % m/s
m_dot_air_cold_per_channel = Area_cross_sect_cold_per_channel * rho_air * air_speed_cold;
% m_dot_air_cold_per_channel_1 = m_dot_air_cold / num_channels;

% Initial conditions - Hot Side (Air not restricted to channel)
inlet_temp_hot = 308.15;   % K
CFM_fan_hot = 59;              % CubicFt/min
volumetric_flow_rate_hot = CFM_fan_hot * ((0.3048^3) / 60);   % m^3/s - conversion factor
m_dot_air_hot = volumetric_flow_rate_hot / rho_air;
fan_area_hot = (pi * 0.04^2 - pi * 0.015^2);                                       % CHANGEME
air_speed_hot = volumetric_flow_rate_hot / fan_area_hot ;         % m/s

% Compute convective coefficient & fin efficiencies
[R_ku_cold, h_cold] = compute_convective_coefficient_cold_NTU(air_speed_cold, area_per_channel, fin_width_cold, Dh_cold_per_channel, m_dot_air_cold_per_channel);
[R_ku_hot, h_hot] = compute_convective_coefficient_hot_without_NTU(air_speed_hot, fin_area_total_hot, fin_width_hot);
overall_fin_eff_hot = compute_fin_efficiency(h_hot, k_fin_hot, fin_thickness_hot, fin_length_hot, num_fins_hot, per_fin_area_hot, fin_area_total_hot);        
overall_fin_eff_cold = compute_fin_efficiency(h_cold, k_fin_cold, fin_thickness_cold, fin_length_cold, num_fins_cold, per_fin_area_cold, fin_area_total_cold);        


%% Print initialization message
fprintf('<strong>***Initialization***\n</strong>');
fprintf('Inlet Air Temperature - Cold Side (T_in_cold): %.3f K \n', inlet_temp_cold);
fprintf('Inlet Air Speed - Cold Side (U_cold): %.1f m/s \n', air_speed_cold);
fprintf('Inlet Air Temperature - Hot Side (T_in_hot): %.3f K \n', inlet_temp_hot);
fprintf('Inlet Air Speed - Hot Side (U_hot): %.1f m/s \n', air_speed_hot);
fprintf('Convective Coefficient Resistance PER CHANNEL (R_ku_c) - Cold Side: %.3f K/W\n', R_ku_cold);
fprintf('Convective Coefficient Resistance (R_ku_h) - Hot Side: %.3f K/W\n', R_ku_hot);
fprintf('Conductive Coefficient Resistance (R_k_hc): %.3f K/W\n\n', R_k_hc);

%% Main Calculation Body

% Initialize data structure for plotting
cooling_power_arr = zeros(J_iters, 1);
heating_power_arr = zeros(J_iters, 1);
power_required_arr = zeros(J_iters, 1);
outlet_temp_cold_arr = zeros(J_iters, 1);
COP_arr = zeros(J_iters, 1);
delta_J_arr = linspace(0, J_max, J_iters);

% Optimal Variables for optimal current
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
    inlet_temp_cold = 308.15;       % Re-initialize inlet_temp
    fprintf('\n<strong>===Iteration %d===\n</strong>', i);
    fprintf('Input Current (J_e): %.2f A \n\n', J_e);
    
    % Accumulated variables for through "num_stages" stages
    cooling_power_total = 0;
    heating_power_total = 0;
    power_required_total = 0;
    
    % Loop based on number of stages
    for j = 1:num_stages
        
        fprintf('<strong>-Stage %d- \n</strong>', j);
        fprintf('Inlet Air Temperature - Cold Side (T_in_cold): %.3f K \n', inlet_temp_cold);
    
        % x = T_h, y = T_c, z = Q_c
        syms x y z
        eqn1 = ((x - y) / R_k_hc) + (overall_fin_eff_hot * (x-inlet_temp_hot) / R_ku_hot) == (num_semi_cond * alpha_seeback * J_e * x) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
        eqn2 = (-(x - y) / R_k_hc) + z == (-num_semi_cond * alpha_seeback * J_e * y) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
        eqn3 = z == (overall_fin_eff_cold * num_channels * (y - inlet_temp_cold) ) / R_ku_cold;

        sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
        T_h_peltier = double(sol.x);
        T_c_peltier = double(sol.y);
        Q_c_peltier = double(sol.z);        % Already factored in cold efficiency and all channels...

        Q_h_peltier = overall_fin_eff_hot * (T_h_peltier - inlet_temp_hot) / R_ku_hot;
        power_conduction_peltier = (T_h_peltier - T_c_peltier) / R_k_hc;
        outlet_temp_cold = inlet_temp_cold + ( (Q_c_peltier / num_channels) / (m_dot_air_cold_per_channel * Cp_air) );
        outlet_temp_hot = inlet_temp_hot + Q_h_peltier/(m_dot_air_hot * Cp_air);
        power_required = num_semi_cond * ((R_e_hc * J_e^2) + (alpha_seeback * J_e * (T_h_peltier - T_c_peltier)) );
        coefficient_performance = -100 * Q_c_peltier / power_required;
        
        % Increment data for variables tracking overall performance
        cooling_power_total = cooling_power_total + Q_c_peltier;
        heating_power_total = heating_power_total + Q_h_peltier;
        power_required_total = power_required_total + power_required;
        
        % Update inlet temperature to outlet temp of previous iteration
        inlet_temp_cold = outlet_temp_cold;
        fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.3f K \n\n', outlet_temp_cold);
        
    end
    
   % Store overall performance from multi-stage analysis
    cooling_power_arr(i) = cooling_power_total;
    heating_power_arr(i) = heating_power_total;
    power_required_arr(i) = power_required_total;
    outlet_temp_cold_arr(i) =  outlet_temp_cold - 273.15;
    
    % Calculate overall COP
    COP_total = -(100 * cooling_power_total) / power_required_total;
    COP_arr(i)= COP_total;
    
    % Find optimal current which gives max cooling
    if -cooling_power_total > -max_cooling_power
        max_cooling_power = cooling_power_total;
        max_heating_power = heating_power_total;
        J_optimal = J_e;
        T_h_optimal = T_h_peltier;
        T_c_optimal = T_c_peltier;
        power_required_optimal = power_required_total;
        COP_optimal = COP_total;
        outlet_temp_cold_optimal = outlet_temp_cold;
        outlet_temp_hot_optimal = outlet_temp_hot;
    end
    
    % Print results..
    fprintf('<strong>Summary \n</strong>');
    fprintf('Overall Power Required (P_e): %.1f W\n', power_required_total);
    fprintf('Coefficient of Performance (COP): %.1f %% \n', COP_total);
    
    % Cold side
    fprintf('---\n');
    fprintf('Cold side Temperature (T_c): %.1f K \n', T_c_peltier);
    fprintf('Overall Cooling Power - Cold Side (Q_c_peltier): %.2f W\n', cooling_power_total);
    fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.1f K\n', outlet_temp_cold);
    
    % Hot side
    fprintf('---\n');
    fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_peltier);
    fprintf('Overall Heating Power - Hot Side (Q_h_peltier): %.2f W\n', heating_power_total);
    fprintf('Outlet Air Temperature - Hot Side (T_out_hot): %.1f K\n\n', outlet_temp_hot);

end

% Display optimal results
% outlet_temp_cold_optimal = inlet_temp_cold + max_cooling_power/(m_dot_air_cold * Cp_air);

fprintf('<strong>===Optimal Current Results===\n</strong>');

fprintf('Input Current (J_e): %.2f A \n', J_optimal);
fprintf('Power Required (P_e): %.1f W\n', power_required_optimal);
fprintf('Coefficient of Performance (COP): %.1f %% \n', COP_optimal);

% Cold side
fprintf('---\n');
fprintf('Cold side Temperature (T_c): %.1f K \n', T_c_optimal);
fprintf('Cooling Power - Cold Side (Q_c_peltier): %.2f W\n', max_cooling_power);
fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.1f K\n', outlet_temp_cold_optimal);

% Hot side
fprintf('---\n');
fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_optimal);
fprintf('Heating Power - Hot Side (Q_h_peltier): %.2f W\n', max_heating_power);
fprintf('Outlet Air Temperature - Hot Side (T_out_hot): %.1f K\n\n', outlet_temp_hot_optimal);


%% Plot final graphs

% Remove data points for high COP for visibility sake
for i = 1:length(delta_J_arr)
    if delta_J_arr(i) < 1.3
        COP_arr(i) = NaN;
    end
end

% % Plot graph of Cooling power against Current
% figure(1)
% plot(delta_J_arr, cooling_power_arr, delta_J_arr, heating_power_arr, delta_J_arr, power_required_arr, delta_J_arr, COP_arr);
% title("Power against Input Current");
% xlabel("Current [A]");
% ylabel("Power [W]");
% legend("Cooling Power (Q_c)", "Heating Power (Q_h)", "Power Input", "COP [%]",  "Location", "NorthEast");
% grid on;
% set(gca,'FontSize',12)

% Plot abs cooling power and power consumption against current

hold on;
figure(1)
plot(delta_J_arr, -cooling_power_arr, delta_J_arr, power_required_arr, delta_J_arr, outlet_temp_cold_arr, delta_J_arr, COP_arr);
title("Input Current Analysis (TEC1-12710) - " + num_stages + " Peltier Modules in Series");
xlabel("Current [A]");
ylabel("Magnitude");
legend("Cooling Power [W]", "Power Consumed [W]", "Outlet Temp [degC]", "COP [%]", "Location", "NorthEast");
grid on;
set(gca,'FontSize',12)

%% Main Functions Used


% Forced convection - Re > 2300 = turbulent flow
% Calculations is per channel...
function [R_ku, h] = compute_convective_coefficient_cold_NTU(air_speed, Area_fin_total, fin_width, Dh, mass_dot)
    
    global kin_visc_air k_air Pr_air Cp_air;
    Re = (air_speed * Dh)/kin_visc_air;
    Nu = 0.023 * Re^(4/5) * Pr_air^(0.3);       % n = 0.3
    NTU = (Area_fin_total * Nu * k_air) / (Dh * mass_dot * Cp_air);
    heat_transfer_eff = 1 - exp(-NTU);                      % R changes with flow in channel...
    R_ku = 1 / (mass_dot * Cp_air * heat_transfer_eff);     % Average convective resistance
%     h = 1 / (R_ku * Area_fin_total);                        % R = 1/hA
    h = Nu * k_air / Dh;
    
%     Re = (air_speed * fin_width)/kin_visc_air;
%     Nu = 0.664 * Re^(0.5) * Pr_air^(1/3);
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



