function [dr_o, ddr, dxi, dw] = se3_dynamics_eul(r, dr, xi, w, fz, tau, d, M, g, J)
%% se3_dynamics: translational-rotational dynamics for 3d body
%  r     [3x1] position            [m]
%  dr    [3x1] velocity            [m/s]
%  w     [3x1] euler angles        [rad]  body to world XYZ
%  w     [3x1] angular velocities  [rad/s] 
%  fz    [1x1] vertical thrust     [N]
%  tau   [3x1] torques             [Nm]
%  d     [3x1] external force      [N]
%  M     [1x1] total weight        [kg]
%  g     [1x1] gravity             [kgm/s^2]
%  J     [3x3] moment of inertia   [g*cm^2] 

  % Euler angles  (body to world XYZ)
  phi = xi(1);
  the = xi(2);
  psi = xi(3);

  % Rotation r_World = Rz*Ry*Rx*r_Body
  R = [cos(psi)*cos(the), cos(psi)*sin(phi)*sin(the) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(the);
       cos(the)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(the), cos(phi)*sin(psi)*sin(the) - cos(psi)*sin(phi);
      -sin(the),          cos(the)*sin(phi),                              cos(phi)*cos(the)];

  dr_o = dr;
  ddr = (1/M) * fz * R(:,3) - [0; 0; g] + d/M;
  [dxi, dw] = so3_eul(xi, w, tau, g, J);
end