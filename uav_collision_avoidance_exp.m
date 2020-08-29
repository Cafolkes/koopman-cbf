
addpath('../uav_sim_ros/codegen/','../uav_sim_ros/codegen/dynamics/','dynamics', 'controllers','koopman_learning','utils')
N = 2;
ts = 1e-3;
Ts = 3e-2;

% Define system and dynamics:
config = quad1_constants;
Kpxy = 2; %0.7
Kpz = 2; %1
KpVxy = 0.7; %0.7
KpVz = 1; %1
KpAtt = 10; %10
KdAtt = 1; %1
KpOmegaz = 2; %2
V_max = 14.8;
V_min = 0.4;
hoverT = 0.5126*V_max; %0.52
Omega_hover = 497.61*ones(4,1);
maxPosErr = 0.25;
M = [Kpxy; Kpz; KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % Backup controller parameters

stop_crit = @(x,x_f) norm(x(:,1)-x_f(:,1)) <= 1e-1 && norm(x(:,2)-x_f(:,2)) <= 1e-1;
legacy_controller = @(x, x_f) pdU(x,x_f,M); % TODO: implement legacy controller.
affine_dynamics = @(x) UAVDynamics(x);                  % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config);     % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                % Processing of state data while simulating

x0_1 = [[1; 1; 1]; eul2quat(zeros(1,3))'; zeros(6,1); Omega_hover];
x0_2 = [[0; 0; 1]; eul2quat(zeros(1,3))'; zeros(6,1); Omega_hover];
xf_1 = [[0; 0; 1]; eul2quat(zeros(1,3))'; zeros(6,1); Omega_hover];
xf_2 = [[1; 1; 1]; eul2quat(zeros(1,3))'; zeros(6,1); Omega_hover];
x0 = [x0_1 x0_2];
xf = [xf_1 xf_2];

%% Run experiment:
[tt,X,U] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, stop_crit, ts, maxPosErr);

%% Plot experiment:
plot_uav_exp(tt,X,U,Ts)

function [tt,X,U] = simulate_sys(x0, xf, sim_dynamics, sim_process, controller, controller_process, stop_criterion, ts, maxPosErr)
    t = 0;
    tt = 0;
    n_agents = size(x0,2);
    x = x0;
    for i = 1 : n_agents
        X{i} = [];
        U{i} = [];
    end
    
    while ~stop_criterion(x,xf)
        for i = 1 : n_agents
            p_d = (xf(1:3,i)-x(1:3,i))*maxPosErr + x(1:3,i);
            x_d = [p_d;xf(4:end,i)];
            u = controller(x(:,i),x_d);
            u = controller_process(u);
            xdot = @(t,x) sim_dynamics(x,u);
            [~, x_tmp] = ode45(xdot,[t,t+ts],x(:,i));
            x(:,i) = x_tmp(end,:)';
            x(:,i) = sim_process(x(:,i), ts);
            X{i} = [X{i};x(:,i)'];
            U{i} = [U{i};u'];
        end
        t = t + ts;
        tt = [tt;t];
    end
end

function plot_uav_exp(tt,X,U,Ts)
    n_agents = size(X,2);
    [x_s, y_s, z_s] = sphere;
    r_s = 0.15;
    for j = 1 : n_agents
       t_int = tt(1):Ts:tt(end);
       X{j} = interp1(tt(1:end-1),X{j},t_int);
    end
    n_data = length(t_int);

    %global Ts am T_exp obs r
    figure('Name','uav simulation','Position', [10 10 900 600])
    for i=1:n_data
        clf
        hold on
        grid on
        for j = 1 : n_agents
            x = X{j}(i,:);
            plotTransforms(x(1:3), x(4:7), 'MeshFilePath', 'multirotor.stl', 'MeshColor','k','FrameSize',0.3)
            plot3(x(1),x(2),x(3),'rO');
            plot3(X{j}(1:i,1),X{j}(1:i,2),X{j}(1:i,3))
            sphere_surface = surf(x(1)+x_s*r_s,x(2)+y_s*r_s,x(3)+z_s*r_s);
            if dist(X{1}(i,1:3)', X{2}(i,1:3)') <= r_s
                 set(sphere_surface,'FaceColor',[1 0 0], ...
                'FaceAlpha',0.3,'FaceLighting','flat','EdgeAlpha',0.3)
            else
                set(sphere_surface,'FaceColor',[0 0 0], ...
                'FaceAlpha',0.1,'FaceLighting','flat','EdgeAlpha',0.3)
            end
            %quiver(x(1),x(2),0.2*x(3)*cos(x(4)),0.2*x(3)*sin(x(4)));
            axis([-0.25 1.25 -0.25 1.25 0 2.5])
            view([0.5 2 1])
        end
        drawnow
    end
end
function d = dist(p1,p2)
    d = (p1-p2)'*(p1-p2);
end