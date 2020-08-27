function [T,X] = collect_data(sim_dynamics, sim_process, controller, controller_process, stop_criterion, initial_condition, n_samples, ts)
    global Ts;
    
    if nargin<8
        ts = 0.02;
    end
    disp('Collecting data...')
    X = {};
    i = 1;
    while i <= n_samples
        x = initial_condition();
        [tt,xx] = simulate_sys(x,sim_dynamics,sim_process,controller, controller_process, stop_criterion, ts);  % Simulate backup trajectory from initial value
        if stop_criterion(tt(end),xx(end,:)) && length(tt) > 2
            xx1 = interp1(tt,xx,0:Ts:tt(end));
            T{i} = tt;
            X{i} = xx1;
            i = i+1;
            fprintf('Sample %i generated\n', i-1);
        end
    end
end

function [tt,xx] = simulate_sys(x, sim_dynamics, sim_process, controller, controller_process, stop_criterion, ts)
    global T_max
    t = 0;
    tt = 0;
    xx = x';
    
    while ~stop_criterion(t,x) && t < T_max
        u = controller(x);
        u = controller_process(u);
        xdot = @(t,x) sim_dynamics(x,u);
        [~, x_tmp] = ode45(xdot,[t,t+ts],x);
        x = x_tmp(end,:)';
        x = sim_process(x, ts);
        t = t + ts;
        tt = [tt;t];
        xx = [xx;x'];
        %disp([x(1:3) x(8:10)])
    end
end