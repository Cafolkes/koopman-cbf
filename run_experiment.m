function [x_rec, u_rec, u0_rec] = run_experiment(x0, system_dynamics, legacy_controller, supervisory_controller)
    %State is defined as x = [X,Y,v,theta], u = [a,r]
    global Ts vm am T_exp

    x = x0;
    x_rec = [];
    u_rec = [];
    u0_rec = [];
   
    for i=1:round(T_exp/Ts)
        u0 = legacy_controller(x);
        
        N = ceil(x(3)/am/Ts);
        u = supervisory_controller(x,u0,N);
    
        [f,g] = system_dynamics(x);
        xdot = f+g*u;
        if x(3)<0
            xdot(3)=max(u(1),-x(3));
        elseif x(3)>vm
            xdot(3)=min(u(1),vm-x(3));
        end
        x = x+xdot*Ts;
        x_rec = [x_rec;x'];
        u_rec = [u_rec;u'];
        u0_rec = [u0_rec;u0'];     
    end
end