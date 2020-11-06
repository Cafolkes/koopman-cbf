function [x_rec, u_rec, u0_rec, comp_t_rec] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller)
    %State is defined as x = [X,Y,v,theta], u = [a,r]
    global Ts am T_exp

    x = x0;
    t = 0;
    x_rec = [];
    u_rec = [];
    u0_rec = [];
    comp_t_rec = [];
   
    for i=1:round(T_exp/Ts)
        u0 = legacy_controller(x);
        
        comp_t0 = posixtime(datetime('now'));
        N = ceil(x(3)/am/Ts);
        u = supervisory_controller(x,u0,N);
        comp_tf = posixtime(datetime('now')) - comp_t0;
        xdot = @(t,x) sim_dynamics(x,u);
        [~, x_tmp] = ode45(xdot,[t,t+Ts],x);
        x = x_tmp(end,:)';
        x = sim_process(x, Ts);

        t = t + Ts;
        x_rec = [x_rec;x'];
        u_rec = [u_rec;u'];
        u0_rec = [u0_rec;u0'];     
        comp_t_rec = [comp_t_rec comp_tf];
    end
end