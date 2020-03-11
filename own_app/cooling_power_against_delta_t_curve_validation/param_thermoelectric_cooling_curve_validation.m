%% Properties for calculation

clear global;
clc;clear;

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global R_e_hc R_k_hc alpha_seeback num_semi_cond


%% Properties



% Table C.9 - Bismuth Telluride Peltier element properties (TEC1-12706)
alpha_s_pos = 2.3 * 10^-4;           % V/degC
alpha_s_neg = -2.1 * 10^-4;          % V/degC
rho_e_pos =  10^-5;                   % Ohm-m
rho_e_neg = rho_e_pos;               % Ohm-m
k_bismuth_pos = 1.7;                 % W/mK
k_bismuth_neg = 1.45;                % W/mK
width_semi_cond = 0.00138;            % m (NOT TOO SURE..)
height_semi_cond = 0.0032;           % m
num_semi_cond = 254;                 % 127 couples

alpha_seeback = alpha_s_pos - alpha_s_neg;
R_e_hc = (height_semi_cond/(width_semi_cond^2)) * (rho_e_pos + rho_e_neg);
R_k_hc = 1 / ( num_semi_cond * (width_semi_cond^2/height_semi_cond) * (k_bismuth_pos + k_bismuth_neg) );

%% Plot cooling power against delta T curve

% T_h = 27 + 273.16;
T_h = 50 + 273.16;

delta_T_max = 80;
iters = 80;
delta_T_arr = linspace(0, delta_T_max, iters);
cooling_power_arr_1_amp = zeros(iters, 1);
cooling_power_arr_2_amp = zeros(iters, 1);
cooling_power_arr_3_amp = zeros(iters, 1);
cooling_power_arr_4_amp = zeros(iters, 1);
cooling_power_arr_5_amp = zeros(iters, 1);
cooling_power_arr_6_amp = zeros(iters, 1);




for i = 1:length(delta_T_arr)
    
    T_c = T_h - delta_T_arr(i);
    
    % 1 A
    J_e = 1;
    cooling_power_arr_1_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
    % 2 A
    J_e = 2;
    cooling_power_arr_2_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
    % 3 A
    J_e = 3;
    cooling_power_arr_3_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
    % 4 A
    J_e = 4;
    cooling_power_arr_4_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
    % 5 A
    J_e = 5;
    cooling_power_arr_5_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
    % 6 A
    J_e = 6;
    cooling_power_arr_6_amp(i) = -(num_semi_cond * alpha_seeback * J_e * T_c) + (delta_T_arr(i) / R_k_hc) + (0.5 * num_semi_cond * R_e_hc * J_e^2);
    
end

plot(delta_T_arr, -cooling_power_arr_1_amp, delta_T_arr, -cooling_power_arr_2_amp, delta_T_arr, -cooling_power_arr_3_amp, delta_T_arr, -cooling_power_arr_4_amp, delta_T_arr, -cooling_power_arr_5_amp, delta_T_arr, -cooling_power_arr_6_amp);
title("Cooling Power against Delta Temp (Hot vs Cold Side)");
xlabel("Delta T [K]");
ylabel("Power [W]");
ylim([1,80]);
legend("1A", "2A", "3A", "4A", "5A", "6A", "Location", "NorthEast");
grid on;

