%% quad1_constants: contains constants for UAV 1
function [con] = quad1_constants(arg)

%   con = struct('m', 2.746, ...                                   % [kg]  everything except motors and props
%                'g', 9.80665, ...                                 % [m/s2]
%                'J_bod', [1e-7*136440 0 0;0 1e-7*142300 0;0 0 1e-7*235160], ...  % [kg m^2] everything except motors and props
%                'D', 0.4766, ...                                  % [m] diagonal length prop-to-prop 
%                'm_z', 0.011, ...                                 % [m] z location of props
%                'm_prop', 0.032, ...                              % [kg]
%                'J_prop', 1e-7*1217/2, ...                             % [kg m^2]
%                'k_f', 2.2234e-07, ...                                % [N/(rad/s)^2]
%                'k_t',1.4179e-05, ...                                % [Nm/(rad/s)^2]
%                'm_mot', 0.108, ...                               % [kg]
%                'J_rot', 1e-7*1217/2, ...                              % [kg m^2]
%                'K_v', 82.2810, ...                             % [(rad/s)/V]
%                'R', 0.074168864848772, ...                                    % [Ohm]
%                'C_g', [0;0;0]);                                  % [m]

con = struct('m', 2.424, ...                                   % [kg]  everything except motors and props
               'g', 9.80665, ...                                 % [m/s2]
               'J_bod', [0.026644 0 0;0 0.029230 0;0 0 0.046516], ...  % [kg m^2] everything except motors and props
               'D', 0.6, ...                                  % [m] diagonal length prop-to-prop 
               'm_z', 0.02, ...                                 % [m] z location of props
               'm_prop', 0.015, ...                              % [kg]
               'J_prop', 3.5e-05, ...                             % [kg m^2]
               'k_f', 2.4e-5, ...                                % [N/(rad/s)^2]
               'k_t',4.5e-7, ...                                % [Nm/(rad/s)^2]
               'm_mot', 0.108, ...                               % [kg]
               'J_rot', 3.5e-05, ...                              % [kg m^2]
               'K_v', 1/0.0136, ...                             % [(rad/s)/V]
               'R', 0.1, ...                                    % [Ohm]
               'C_g', [0;0;0]);                                  % [m]