%% omega_to_torque: body torques and force from motor angular velocities
function [F_z, tau] = body_forces(Omega, D, k_f, k_t, Roc)

  l = D/(2*sqrt(2));

  r1 = [l  -l 0]';
  r2 = [-l -l 0]';
  r3 = [-l  l 0]';
  r4 = [l   l 0]';

  F1 = [0 0 k_f*Omega(1)^2]';
  F2 = [0 0 k_f*Omega(2)^2]';
  F3 = [0 0 k_f*Omega(3)^2]';
  F4 = [0 0 k_f*Omega(4)^2]';

  T1 = -[0 0 k_t*Omega(1)^2]';
  T2 =  [0 0 k_t*Omega(2)^2]';
  T3 = -[0 0 k_t*Omega(3)^2]';
  T4 =  [0 0 k_t*Omega(4)^2]';

  F_tot = F1 + F2 + F3 + F4;
  T_tot = cross(r1, F1) + cross(r2, F2) + cross(r3, F3) + cross(r4, F4) ...
            + T1 + T2 + T3 + T4;

  if nargin >= 5
    T_tot = T_tot - cross(Roc, F_tot);  % compensate for center of gravity
  end

  F_z = F_tot(3);
  tau = T_tot;
end