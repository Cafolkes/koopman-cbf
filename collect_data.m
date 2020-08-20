function X = collect_data(system_dynamics, controller, stop_criterion, n_samples)
    %State is defined as x = [X,Y,v,theta], u = [a,r]
    global Ts x_bdry;

    X = {};

    for i=1:n_samples
        x = x_bdry(:,1)+(x_bdry(:,2)-x_bdry(:,1)).*rand(4,1);  % Sample random initial value of x inside x_bdry
        [tt,xx] = generate_backup_traj(x,system_dynamics,controller,stop_criterion);  % Simulate backup trajectory from intial value
        xx1 = interp1(tt,xx,0:Ts:tt(end));
        X{i} = xx1;
    end
end

function ja = dubin_f_x(x,controller)
    h = 1e-6;
    dudt = zeros(2,4);
    dudt(:,1) = (controller(x+[h;0;0;0])-controller(x-[h;0;0;0]))/2/h;
    dudt(:,2) = (controller(x+[0;h;0;0])-controller(x-[0;h;0;0]))/2/h;
    dudt(:,3) = (controller(x+[0;0;h;0])-controller(x-[0;0;h;0]))/2/h;
    dudt(:,4) = (controller(x+[0;0;0;h])-controller(x-[0;0;0;h]))/2/h;
    ja = [0 0 cos(x(4)) -x(3)*sin(x(4));...
        0 0 sin(x(4)) x(3)*cos(x(4));...
        dudt];
end

function [tt,xx] = generate_backup_traj(x,system_dynamics, controller,stop_criterion,ts)
    global am vm

    if nargin<5
        ts = 0.02;
    end
    t = 0;
    tt = 0;
    xx = x';

    while ~stop_criterion(t,x)

        u = controller(x);
        [f,g] = system_dynamics(x);
        xdot = f+g*u;
        if x(3)<-vm
            xdot(3)=max(u(1),-x(3)-vm);
        elseif x(3)>vm
            xdot(3)=min(u(1),vm-x(3));
        end

        ja = dubin_f_x(x,controller);
        x = x + xdot*ts;
        if abs(x(3))<ts*am/2
            x(3)=0;
        end
        t = t + ts;
        tt = [tt;t];
        xx = [xx;x'];

    end
end