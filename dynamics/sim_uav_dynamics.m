function xdot = sim_uav_dynamics(x,u,config,use_quat,add_ground_effect)
    if use_quat == true
        r = x(1:3);
        q = x(4:7);
        v = x(8:10);
        w = x(11:13);
        Omega = x(14:17);
    else
        r = x(1:3);
        q = x(4:6);
        v = x(7:9);
        w = x(10:12);
        Omega = x(13:16);
    end
    if add_ground_effect == true
        d = ground_effect_force(r,v,q,w,Omega,config.R,sqrt(2*(config.D/2)^2),config.D,config.k_f,config.k_t,use_quat);
    else
        d = zeros(3,1);
    end
    [dr, ddr, dxi_q, dw, dOmega] = full_eul(r,v,q,w,Omega,u,d,config.m,config.g,config.J_bod,config.D,config.J_rot, config.J_prop, config.K_v, config.R, config.k_f, config.k_t,use_quat);
    xdot = [dr;dxi_q;ddr;dw;dOmega];
end

function f_g = ground_effect_force(r,v,q,w,Omega,R,d,D,k_f,k_t,use_quat,K_b,noise)
  if nargin < 12
      K_b = 1;
  end
  if nargin < 13
      noise = 0.02;
  end
  z = r(3);
  den = 1-(R/(4*z))^2 - R^2*(z/sqrt((d^2+4*z^2)^3)) - (R^2/2)*(z/sqrt((2*d^2+4*z^2)^3)) - 2*R^2*(z/sqrt((D^2+4*z^2)^3))*K_b;
  ge_force_factor = 1/den;
  ge_force_factor = ge_force_factor + max(0,1-z)*normrnd(0,noise);
  
  F1 = [0 0 k_f*Omega(1)^2*ge_force_factor]';
  F2 = [0 0 k_f*Omega(2)^2*ge_force_factor]';
  F3 = [0 0 k_f*Omega(3)^2*ge_force_factor]';
  F4 = [0 0 k_f*Omega(4)^2*ge_force_factor]';
  F_tot = F1 + F2 + F3 + F4;
  
  if use_quat == true
      rotm = quat2romt(q);
  else
      %rotm = eul2rotm(q','XYZ');
        phi = q(1);
        the = q(2);
        psi = q(3);
  
        rotm = [cos(psi)*cos(the), cos(psi)*sin(phi)*sin(the) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(the);
              cos(the)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(the), cos(phi)*sin(psi)*sin(the) - cos(psi)*sin(phi);
              -sin(the),          cos(the)*sin(phi),                              cos(phi)*cos(the)];
  end
  
  f_g = rotm*F_tot-rotm*[0;0;body_forces(Omega, D, k_f, k_t)];
end