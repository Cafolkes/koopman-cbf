function [dxi, dw] = se3_dynamics_eul(xi, w, tau, g, J)
%% se3_dynamics: rotational dynamics for 3d body
%  xi    [3x1] euler angles        [rad]  body to world XYZ
%  r     [3x1] rates               [rad/s]
%  tau   [3x1] torques             [Nm]
%  g     [1x1] gravity             [kgm/s^2]
%  J     [3x3] moment of inertia   [kg*m^2] 

  % Euler angles  (body to world XYZ)
  phi = xi(1);
  the = xi(2);
  
  secant_fn = @(x) 1/cos(x); 

  T = [1 sin(phi)*tan(the)  cos(phi)*tan(the);
       0 cos(phi)          -sin(phi);
       0 secant_fn(the)*sin(phi)  cos(phi)*secant_fn(the)];

  dxi = T * w;
  dw = inv(J)*(tau - cross(w, J*w));
end