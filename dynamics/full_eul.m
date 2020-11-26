function [dr, ddr, dxi_q, dw, dOmega] = full_eul(r, dr, xi_q, w, Omega, V, d, ...
                                                      M, g, J, D, J_rot, J_prop, ...
                                                      K_v, R, k_f, k_t, ...
                                                      use_quat)
%% full_dynamics_quat: se3 + motor dynamics
%  r      [3x1] position            [m]
%  dr     [3x1] velocity            [m/s]
%  xi_q   [4x1]/[3x1] Unit quaternion or Euler angles
%  w      [3x1] angular velocities  [rad/s] 
%  Omega  [4x1] motor ang velocity  [rad/s]
%  V      [4x1] input voltage       [V]
%  d      [3x1] external force      [N]
% 
%  M      [1x1] total weight        [kg]
%  g      [1x1] gravity             [kgm/s^2]
%  J      [3x3] moment of inertia   [g*cm^2] 
%  D      [1x1] frame diameter      [m]
%  J_rot  [1x1] Inertia rotors      [g cm^2]
%  J_prop [1x1] Inertia propellers  [g cm^2]
%
%  K_v     Motor constant           [rpm/V]
%  R       Motor resistance         [Ohm]
%  k_t     Torque constant          [Nm/(rad/s^2)]
%  k_f     Force constant           [Nm/(rad/s^2)]
% 
% use_quat Use quaternions          default:false

  if nargin < 18
    use_quat = false
  end

  % Motor dynamics
  dOmega = [motor(Omega(1), V(1), J_rot + J_prop, K_v, R, k_t);
            motor(Omega(2), V(2), J_rot + J_prop, K_v, R, k_t);
            motor(Omega(3), V(3), J_rot + J_prop, K_v, R, k_t);
            motor(Omega(4), V(4), J_rot + J_prop, K_v, R, k_t)];


  % Motor velocities to body forces
  [F_z, tau] = body_forces(Omega, D, k_f, k_t);

  % Body dynamics
  if use_quat
    [dr, ddr, dxi_q, dw] = se3_quat(r, dr, xi_q, w, F_z, tau, d, M, g, J);
  else
    [dr, ddr, dxi_q, dw] = se3_eul(r, dr, xi_q, w, F_z, tau, d, M, g, J);
  end