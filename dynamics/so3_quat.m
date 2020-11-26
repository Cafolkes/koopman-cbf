function [dq, dw] = so3_quat(q, w, tau, g, J)
%% se3_dynamics: rotational dynamics for 3d body
%  q     [4x1] quaternion          body to world XYZ
%  r     [3x1] rates               [rad/s]
%  tau   [3x1] torques             [Nm]
%  g     [1x1] gravity             [kgm/s^2]
%  J     [3x3] moment of inertia   [kg*m^2] 

  dq = quatmultiply(q', [0 w'])'/2;                 % w velocity of frame w.r.t. world
  dw = inv(J)*(tau - cross(w, J*w));
end