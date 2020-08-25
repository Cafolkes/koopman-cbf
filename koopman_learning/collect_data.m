function X = collect_data(sim_dynamics, sim_process, controller, stop_criterion, n_samples)
    %State is defined as x = [X,Y,v,theta], u = [a,r]
    global Ts x_bdry;

    X = {};

    for i=1:n_samples
        x = x_bdry(:,1)+(x_bdry(:,2)-x_bdry(:,1)).*rand(4,1);  % Sample random initial value of x inside x_bdry
        [tt,xx] = simulate_sys(x,sim_dynamics,sim_process,controller,stop_criterion);  % Simulate backup trajectory from intial value
        xx1 = interp1(tt,xx,0:Ts:tt(end));
        X{i} = xx1;
    end
end

function [tt,xx] = simulate_sys(x, sim_dynamics, sim_process, controller, stop_criterion, ts)
    %global vm am Ts    
    if nargin<6
        ts = 0.02;
    end
    
    t = 0;
    tt = 0;
    xx = x';

    while ~stop_criterion(t,x)
    %for i = 1 : ceil(vm/am/Ts)
        u = controller(x);
        xdot = @(t,x) sim_dynamics(x,u);
        [~, x_tmp] = ode45(xdot,[t,t+ts],x);
        x = x_tmp(end,:)';
        x = sim_process(x, ts);
        t = t + ts;
        tt = [tt;t];
        xx = [xx;x'];
    end
end