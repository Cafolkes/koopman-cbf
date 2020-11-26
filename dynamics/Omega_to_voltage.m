%% Omega_to_voltage: function description
function [V] = Omega_to_voltage(Omega, K_v, k_t, R)
  V = Omega/K_v + (k_t*Omega.^2)*K_v*R;
end