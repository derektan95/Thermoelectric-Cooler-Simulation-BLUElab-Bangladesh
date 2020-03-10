%% Load essential parameters
% Implement struct data structure in the future!

% warning('off','all');           % Turn off all warnings
run("param_thermoelectric_cooling_dimension_analysis.m");

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global Area_cross_sect_cold_per_channel Dh_cold_per_channel num_channels area_per_channel
global R_e_hc R_k_hc alpha_seeback num_semi_cond
global fin_width_cold fin_length_cold fin_thickness_cold sink_height_cold num_fins_cold k_fin_cold per_fin_area_cold base_area_cold fin_area_total_cold
global fin_width_hot fin_length_hot fin_thickness_hot sink_height_hot num_fins_hot k_fin_hot per_fin_area_hot base_area_hot fin_area_total_hot 

%% Define simulation parameters (CHANGME)

% % To iterate fin Width
% fin_width = 0.045;          % Original short fin width
% fin_width_iters = 100;
% fin_width_max = 0.45;       % 10 times the length
% 
% % To iterate fin Length
% fin_length = 0.021;
% fin_length_iters = 100;
% fin_length_max = 0.21;       % 10 times the length
%
% % To iterate sink height
sink_height_extra_channels = 0;
sink_height_extra_channels_iters = 100;
sink_height_extra_channels_max = 100;       % 10 times the length


% General parameters
J_e = 0.93;   

% Initial conditions - Cold Side (Air restricted to channel)
inlet_temp_cold = 308.15;   % K
air_speed_cold = 2.2;      % m/s
m_dot_air_cold_per_channel = Area_cross_sect_cold_per_channel * rho_air * air_speed_cold;

% Initial conditions - Hot Side (Air not restricted to channel)
inlet_temp_hot = 308.15;   % K
CFM_fan_hot = 48;              % CubicFt/min
volumetric_flow_rate_hot = CFM_fan_hot * ((0.3048^3) / 60);   % m^3/s - conversion factor
m_dot_air_hot = volumetric_flow_rate_hot / rho_air;
fan_area = pi * 0.04^2;
air_speed_hot = volumetric_flow_rate_hot / fan_area ;         % m/s


%% Print initialization message
fprintf('<strong>***Initialization***\n</strong>');
fprintf('Inlet Air Temperature - Cold Side (T_in_cold): %.3f K \n', inlet_temp_cold);
fprintf('Inlet Air Speed - Cold Side (U_cold): %.1f m/s \n', air_speed_cold);
fprintf('Inlet Air Temperature - Hot Side (T_in_hot): %.3f K \n', inlet_temp_hot);
fprintf('Inlet Air Speed - Hot Side (U_hot): %.1f m/s \n', air_speed_hot);
% fprintf('Convective Coefficient Resistance PER CHANNEL (R_ku_c) - Cold Side: %.3f K/W\n', R_ku_cold);
% fprintf('Convective Coefficient Resistance (R_ku_h) - Hot Side: %.3f K/W\n', R_ku_hot);
% fprintf('Conductive Coefficient Resistance (R_k_hc): %.3f K/W\n\n', R_k_hc);

%% Main Calculation Body

% delta_fin_width_arr = linspace(fin_width, fin_width_max, fin_width_iters);
% delta_fin_length_arr = linspace(fin_length, fin_length_max, fin_length_iters);
delta_sink_height_extra_channels_arr = linspace(sink_height_extra_channels, sink_height_extra_channels_max, sink_height_extra_channels_iters);

cooling_power_arr = zeros(sink_height_extra_channels_iters, 1);
power_required_arr = zeros(sink_height_extra_channels_iters, 1);
J_optimal = 0;
max_cooling_power = 0;
T_h_optimal = 0;
T_c_optimal = 0;
power_required_optimal = 0;
COP_optimal = 0;
outlet_temp_cold_optimal = 0;
outlet_temp_hot_optimal = 0;

for i = 1:length(delta_sink_height_extra_channels_arr)
    
    % Redefine fin conditions - Cold Side 
    
    % Param that stay constant
    width_btwn_fins = 0.0039;
    k_fin_cold = 237;                % Conduction Coeff - Aluminum [W/mK]
    fin_thickness_cold = 0.001;          % CHANGEME
    
    fin_width_cold = 0.045;           % length parallel to flow [m]
%     fin_width_cold = delta_fin_width_arr(i);  

    fin_length_cold = 0.021;
%     fin_length_cold = delta_fin_length_arr(i);  

    % Don't increment on 1st iteration
    if i > 0
        sink_height_cold = sink_height_cold + width_btwn_fins + fin_thickness_cold;           
        num_fins_cold = num_fins_cold + 1;    
    end
    delta_sink_height_extra_channels_arr(i) = sink_height_cold;
    
    
    per_fin_area_cold = 2 * fin_width_cold * fin_length_cold;
    base_area_cold = (fin_width_cold * sink_height_cold) - (num_fins_cold * fin_width_cold * fin_thickness_cold);  
    fin_area_total_cold = ( (num_fins_cold-1) * per_fin_area_cold) + base_area_cold;
    
    % Calculate flow channel
    height = fin_length_cold;
    num_channels = num_fins_cold - 1;
%     width_btwn_fins = (sink_height_cold - (num_fins_cold * fin_thickness_cold) )/num_channels;
    Area_cross_sect_cold_per_channel = (height * width_btwn_fins);
    perimeter_per_channel = (2 * height) + (2 * width_btwn_fins);
    Dh_cold_per_channel = 4*Area_cross_sect_cold_per_channel/perimeter_per_channel; 

    area_per_channel = fin_area_total_cold / num_channels;
 
    % Compute convective coefficient & fin efficiencies
    [R_ku_cold, h_cold] = compute_convective_coefficient_cold_NTU(air_speed_cold, area_per_channel, fin_width_cold, Dh_cold_per_channel, m_dot_air_cold_per_channel);
    [R_ku_hot, h_hot] = compute_convective_coefficient_hot_without_NTU(air_speed_hot, fin_area_total_hot, fin_width_hot);
    overall_fin_eff_hot = compute_fin_efficiency(h_hot, k_fin_hot, fin_thickness_hot, fin_length_hot, num_fins_hot, per_fin_area_hot, fin_area_total_hot);        
    overall_fin_eff_cold = compute_fin_efficiency(h_cold, k_fin_cold, fin_thickness_cold, fin_length_cold, num_fins_cold, per_fin_area_cold, fin_area_total_cold);        

    
    
    % x = T_h, y = T_c, z = Q_c
    syms x y z
    eqn1 = ((x - y) / R_k_hc) + (overall_fin_eff_hot * (x-inlet_temp_hot) / R_ku_hot) == (num_semi_cond * alpha_seeback * J_e * x) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
    eqn2 = (-(x - y) / R_k_hc) + z == (-num_semi_cond * alpha_seeback * J_e * y) + (0.5 * num_semi_cond * R_e_hc * J_e^2); 
    eqn3 = z == (overall_fin_eff_cold * num_channels * (y - inlet_temp_cold) ) / R_ku_cold;

    sol = solve([eqn1, eqn2, eqn3], [x, y, z]);
    T_h_peltier = double(sol.x);
    T_c_peltier = double(sol.y);
    Q_c_peltier = double(sol.z);
    
    Q_h_peltier = (T_h_peltier - inlet_temp_hot) / R_ku_hot;
    power_conduction_peltier = (T_h_peltier - T_c_peltier) / R_k_hc;
    outlet_temp_cold = inlet_temp_cold + ( (Q_c_peltier / num_channels) / (m_dot_air_cold_per_channel * Cp_air) );
    outlet_temp_hot = inlet_temp_hot + Q_h_peltier/(m_dot_air_hot * Cp_air);
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
    fprintf('Input Current (J_e): %.2f A \n', J_e);
    fprintf('Power Required (P_e): %.1f W\n', power_required);
    fprintf('Coefficient of Performance (COP): %.1f %% \n', coefficient_performance);
    
    % Cold side
    fprintf('---\n');
    fprintf('Cold side Temperature (T_c): %.1f K \n', T_c_peltier);
    fprintf('Cooling Power - Cold Side (Q_c_peltier): %.2f W\n', Q_c_peltier);
    fprintf('Outlet Air Temperature - Cold Side (T_out_cold): %.1f K\n', outlet_temp_cold);
    
    % Hot side
    fprintf('---\n');
    fprintf('Hot side Temperature (T_h): %.1f K \n', T_h_peltier);
    fprintf('Heating Power - Hot Side (Q_h_peltier): %.2f W\n', Q_h_peltier);
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

% Plot graph of Cooling power against Current
% figure(1)
% plot(delta_fin_width_arr, cooling_power_arr);
% title("Cooling Power against Fin Width");
% xlabel("Fin Width [m]");
% ylabel("Cooling Power [W]");
% grid on;

% figure(1)
% plot(delta_fin_length_arr, cooling_power_arr);
% title("Cooling Power against Fin Length");
% xlabel("Fin Length [m]");
% ylabel("Cooling Power [W]");
% grid on;

figure(1)
plot(delta_sink_height_extra_channels_arr, cooling_power_arr);
title("Cooling Power against Sink Height");
xlabel("Sink Height [m]");
ylabel("Cooling Power [W]");
grid on;


% Plot abs cooling power and power consumption against current
% hold on;
% figure(2)
% plot(delta_fin_length_arr, -cooling_power_arr, delta_fin_length_arr, power_required_arr);
% title("Absolute Cooling Power and Power Consumption against Fin Length");
% xlabel("Fin Length [m]");
% ylabel("Power [W]");
% legend("Cooling Power", "Power Consumed", "Location", "NorthEast");
% grid on;

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



