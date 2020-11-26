function [dr_o, ddr, dq, dw] = se3_quat(r, dr, q, w, fz, tau, d, M, g, J)
%% se3_dynamics: translational-rotational dynamics for 3d body
%  r     [3x1] position            [m]
%  dr    [3x1] velocity            [m/s]
%  q     [4x1] Unit quaternion
%  w     [3x1] angular velocities  [rad/s] 
%  fz    [1x1] vertical thrust     [N]
%  tau   [3x1] torques             [Nm]
%  d     [3x1] external force      [N]
%  M     [1x1] total weight        [kg]
%  g     [1x1] gravity             [kgm/s^2]
%  J     [3x3] moment of inertia   [g*cm^2] 

  dr_o = dr;
  ddr = (fz/M) * my_quatrotate(q', [0 0 1])' - [0; 0; g] + d/M;  % q rotation body -> world
  [dq, dw] = so3_quat(q, w, tau, g, J);
end

function [r] = my_quatrotate(q, w)
  % inverse of standard Matlab quatrotate
  r = quatmultiply(quatmultiply(q, [0 w]), quatinv(q));
  r = r(2:4);
end