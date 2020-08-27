function xdot = sim_uav_dynamics(x,u,config)
    [dr, ddr, dxi_q, dw, dOmega] = full_eul(x(1:3),x(4:6),x(7:10),x(11:13),x(14:17),u,zeros(3,1),config.m,config.g,config.J_bod,config.D,config.J_rot, config.J_prop, config.K_v, config.R, config.k_f, config.k_t,true);
    xdot = [dr;dxi_q;ddr;dw;dOmega];
end