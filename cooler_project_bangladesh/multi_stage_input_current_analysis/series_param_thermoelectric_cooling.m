%% Properties for calculation

clear global;
clc;clear;

%% Declare variables as global for use in other scripts (bad practice)
global kin_visc_air Cp_air k_air alpha_air Pr_air rho_air 
global Area_cross_sect_cold_per_channel Dh_cold_per_channel num_channels area_per_channel
global R_e_hc R_k_hc alpha_seeback num_semi_cond
global fin_width_cold fin_length_cold fin_thickness_cold sink_height_cold num_fins_cold k_fin_cold per_fin_area_cold base_area_cold fin_area_total_cold
global fin_width_hot fin_length_hot fin_thickness_hot sink_height_hot num_fins_hot k_fin_hot per_fin_area_hot base_area_hot fin_area_total_hot 

%% Properties

% Table C.22 (300K) - Air properties
kin_visc_air = 1.566 * 10^-5;        % Calculate Reynolds (m^2/s)
Cp_air = 1005;                       % heat capacity (J/kgK)
k_air = 0.0267;                      % To calculate Conv Coeff (W/mK)
alpha_air = 2.257 * 10^-5;           % To calculate Le (m^2/s)
Pr_air = 0.69;                       % Prandt Num ()
rho_air = 1.177;                     % Density (Kg/m^3)


% % Table C.9 - Bismuth Telluride Peltier element properties (TEC1-12706) -
% % match Spec Sheet...
% alpha_s_pos = 1.25 * 10^-4;           % V/degC
% alpha_s_neg = -1.25 * 10^-4;          % V/degC
% rho_e_pos = 0.28 * 10^-5;                   % Ohm-m
% rho_e_neg = rho_e_pos;               % Ohm-m
% k_bismuth_pos = 2.0;                 % W/mK
% k_bismuth_neg = 2.0;                % W/mK
% width_semi_cond = 0.0011;            % m 
% height_semi_cond = 0.0022;           % m
% num_semi_cond = 220;                 % 127 couples

% % OLD
% alpha_s_pos = 1.05 * 10^-4;           % V/degC
% alpha_s_neg = -1.05 * 10^-4;          % V/degC
% rho_e_pos = 0.5 * 10^-5;                   % Ohm-m
% rho_e_neg = rho_e_pos;               % Ohm-m
% k_bismuth_pos = 0.8;                 % W/mK
% k_bismuth_neg = 0.8;                % W/mK
% width_semi_cond = 0.0018;            % m (NOT TOO SURE..)
% height_semi_cond = 0.0027;           % m
% num_semi_cond = 254;                 % 127 couples

% % USED FOR MORE POWER BUT MORE INPUT CURRENT...
% % Bismuth Telluride Peltier element properties (TEC1-12710)
alpha_s_pos = 1.15 * 10^-4;           % V/degC
alpha_s_neg = -1.1 * 10^-4;          % V/degC
rho_e_pos = 0.28 * 10^-5;                   % Ohm-m
rho_e_neg = rho_e_pos;               % Ohm-m
k_bismuth_pos = 1.7;                 % W/mK
k_bismuth_neg = 1.65;                % W/mK
width_semi_cond = 0.00135;            % m 
height_semi_cond = 0.0019;           % m
num_semi_cond = 252;                 % 127 couples

alpha_seeback = alpha_s_pos - alpha_s_neg;
R_e_hc = (height_semi_cond/(width_semi_cond^2)) * (rho_e_pos + rho_e_neg);
R_k_hc = 1 / ( num_semi_cond * (width_semi_cond^2/height_semi_cond) * (k_bismuth_pos + k_bismuth_neg) );


%% Fin conditions - Cold Side      

% % SMALL COLD FIN (Fin 1)
% fin_width_cold = 0.045;           % length parallel to flow [m]
% fin_length_cold = 0.021;             % CHANGEME
% fin_thickness_cold = 0.001;          % CHANGEME
% sink_height_cold = 0.04;            % CHANGEME
% num_fins_cold = 9;               % CHANGEME
% k_fin_cold = 237;                % Conduction Coeff - Aluminum [W/mK]

% % LARGE COLD FIN (Same as on hot side) - (Fin 2)
% fin_width_cold = 0.12;                % length parallel to flow [m]
% fin_length_cold = 0.027;                              % CHANGEME
% fin_thickness_cold = 0.001;                           % CHANGEME
% sink_height_cold = 0.08;                             % CHANGEME
% num_fins_cold = 12;                                % CHANGEME
% k_fin_cold = 237;                                 % Conduction Coeff - Aluminum [W/mK]

% % Moderate Length FIN (Fin 3)
% fin_width_cold = 0.13;           % length parallel to flow [m]
% fin_length_cold = 0.03125;             % CHANGEME
% fin_thickness_cold = 0.0007;          % CHANGEME
% sink_height_cold = 0.0693;            % CHANGEME
% num_fins_cold = 27;               % CHANGEME
% k_fin_cold = 237;                % Conduction Coeff - Aluminum [W/mK]

% % Longer Length FIN (Fin 4)
fin_width_cold = 0.15;           % length parallel to flow [m]
fin_length_cold = 0.03125;             % CHANGEME
fin_thickness_cold = 0.0007;          % CHANGEME
sink_height_cold = 0.0693;            % CHANGEME
num_fins_cold = 27;               % CHANGEME
k_fin_cold = 237;                % Conduction Coeff - Aluminum [W/mK]

% % Shorter Length FIN (Fin 5) - Esp for 2 stage analysis
% fin_width_cold = 0.10;           % length parallel to flow [m]
% fin_length_cold = 0.03125;             % CHANGEME
% fin_thickness_cold = 0.0007;          % CHANGEME
% sink_height_cold = 0.0693;            % CHANGEME
% num_fins_cold = 27;               % CHANGEME
% k_fin_cold = 237;                % Conduction Coeff - Aluminum [W/mK]

per_fin_area_cold = 2 * fin_width_cold * fin_length_cold;
base_area_cold = (fin_width_cold * sink_height_cold) - (num_fins_cold * fin_width_cold * fin_thickness_cold);  
fin_area_total_cold = ( (num_fins_cold-1) * per_fin_area_cold) + base_area_cold;
% % fin_area_total_cold = 0.05;      % Given by prof's example [m^2]

  


%% Fin conditions - Hot Side     

% % Defualt hot fins
% fin_width_hot = 0.12;                % length parallel to flow [m]
% fin_length_hot = 0.027;                              % CHANGEME
% fin_thickness_hot = 0.001;                           % CHANGEME
% sink_height_hot = 0.08;                             % CHANGEME
% num_fins_hot = 12;                                % CHANGEME
% k_fin_hot = 237;                                 % Conduction Coeff - Aluminum [W/mK]
    
% % Longer Length FIN (Fin 4)
% fin_width_hot = 0.15;           % length parallel to flow [m]
% fin_length_hot = 0.03125;             % CHANGEME
% fin_thickness_hot = 0.0007;          % CHANGEME
% sink_height_hot = 0.0693;            % CHANGEME
% num_fins_hot = 27;               % CHANGEME
% k_fin_hot = 237;                % Conduction Coeff - Aluminum [W/mK]

% % Fin 5
% fin_width_hot = 0.10;           % length parallel to flow [m]
% fin_length_hot = 0.03125;             % CHANGEME
% fin_thickness_hot = 0.0007;          % CHANGEME
% sink_height_hot = 0.0693;            % CHANGEME
% num_fins_hot = 27;               % CHANGEME
% k_fin_hot = 237;                % Conduction Coeff - Aluminum [W/mK]

% % Cut Fin 4 (150mm) into half
fin_width_hot = 0.075;           % length parallel to flow [m]
fin_length_hot = 0.03125;             % CHANGEME
fin_thickness_hot = 0.0007;          % CHANGEME
sink_height_hot = 0.0693;            % CHANGEME
num_fins_hot = 27;               % CHANGEME
k_fin_hot = 237;                % Conduction Coeff - Aluminum [W/mK]


per_fin_area_hot = 2 * fin_width_hot * fin_length_hot;                              
base_area_hot = (fin_width_hot * sink_height_hot) - (num_fins_hot * fin_width_hot * fin_thickness_hot); 
fin_area_total_hot = ( (num_fins_hot-1) * per_fin_area_hot) + base_area_hot;
% fin_area_total_hot = fin_area_total_cold*4;      % Given by prof's example [m^2]  


%% Dimensions of flow channel (between sheets)
% According to size of heat sink on cold side... (account for area covered
% by fins)

% width = sink_height_cold;
height = fin_length_cold;
num_channels = num_fins_cold - 1;
width_btwn_fins = (sink_height_cold - (num_fins_cold * fin_thickness_cold) )/num_channels;
Area_cross_sect_cold_per_channel = (height * width_btwn_fins);
perimeter_per_channel = (2 * height) + (2 * width_btwn_fins);
Dh_cold_per_channel = 4*Area_cross_sect_cold_per_channel/perimeter_per_channel; 

area_per_channel = fin_area_total_cold / num_channels;

