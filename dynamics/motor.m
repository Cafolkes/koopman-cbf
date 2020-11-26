function [dOmega] = motor_dynamics(Omega, V, J, K_v, R, k_t)
%% motor_dynamics: differential equation for motor/propeller
%  Omega [1x1]  Angular velocity       [rad/s]
%  V     [1x1]  Input voltage          [V]
%  J     [1x1]  Inertia rotating part  [g cm^2]
%  K_v   [1x1]  Motor constant         [(rad/s)/V]
%  R     [1x1]  Motor resistance       [Ohm]
%  k_t   [1x1]  Torque constant        [Nm/(rad/s^2)]
  dOmega = ((1/(K_v*R)) * (V - Omega/(K_v)) - k_t   * Omega^2)/J;
end