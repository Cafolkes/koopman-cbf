function xdot = sim_uav_dynamics(x,u,config,use_quat)
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
    [dr, ddr, dxi_q, dw, dOmega] = full_eul(r,v,q,w,Omega,u,zeros(3,1),config.m,config.g,config.J_bod,config.D,config.J_rot, config.J_prop, config.K_v, config.R, config.k_f, config.k_t,use_quat);
    xdot = [dr;dxi_q;ddr;dw;dOmega];
end