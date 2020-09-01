
addpath('../uav_sim_ros/codegen/','../uav_sim_ros/codegen/dynamics/','dynamics', 'controllers','koopman_learning','utils','utils/qpOASES-3.1.0/interfaces/matlab/')
file_name = 'data/uav_collision_avoidance_eul.mat';               % File to save data matrices
N = 3;
ts = 1e-3;
Ts = 3e-2;

% Define system and dynamics:
config = quad1_constants;
Kpxy = 0.7; %0.7
Kpz = 1; %1
KpVxy = 0.7; %0.7
KpVz = 1; %1
KpAtt = 10; %10
KdAtt = 1; %1
KpOmegaz = 2; %2
V_max = 14.8;
V_min = 0.5;
u_lim = [V_min*ones(4,1) V_max*ones(4,1)];
hoverT = 0.5126*V_max; %0.52
Omega_hover = 497.61*ones(4,1);
maxPosErr = 0.25;
M = [Kpxy; Kpz; KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % PD controller parameters
M_backup = [KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % Backup controller parameters

stop_crit = @(x,x_f) all(x(3,:) <= z_land*ones(1,3));
legacy_controller = @(x, x_f) pdU_eul(x,x_f,M); 
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config,false);     % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                % Processing of state data while simulating


x0_1 = [[sin(2*pi/3); cos(2*pi/3); 2]; zeros(9,1); Omega_hover];
x0_2 = [[sin(4*pi/3); cos(4*pi/3); 2]; zeros(9,1); Omega_hover];
x0_3 = [[sin(0); cos(0); 2]; zeros(9,1); Omega_hover];

xf_1 = [[0; 0; 0]; zeros(9,1); Omega_hover];
xf_2 = [[0; 0; 0]; zeros(9,1); Omega_hover];
xf_3 = [[0; 0; 0]; zeros(9,1); Omega_hover];
x0 = [x0_1 x0_2 x0_3];
xf = [xf_1 xf_2 xf_3];
z_land = 0.1;

% Define Koopman supervisory controller:
%koopman_file = 'data/uav_learned_koopman_eul.mat';         % File containing learned Koopman model
r_margin = 0.15;                                    % Minimum distance between robot center points                
alpha = 0.5;                                          % CBF strengthening term

%load(koopman_file)
func_dict = @(x) uav_D_full(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17));
options = qpOASES_options('printLevel',0);              % Solver options for supervisory controller
affine_dynamics = @(x) UAVDynamics_eul(x);                  % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                        % State is defined as x = [p,q,v,w,Omega], u = [V1,V2,V3,V4]
barrier_func = @(x1,x2) collision_avoidance_3d_vec(x1,x2,r_margin);             % Barrier function
backup_controller = @(x) backupU_eul(x,M_backup);
backup_dynamics = @(x) cl_dynamics(x, affine_dynamics, backup_controller); 

supervisory_controller = @(x,u0,agent_ind) koopman_qp_cbf_multiagent_vec_mod(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func, alpha, N, func_dict, K_pows, C, options, u_lim,16,4);

%% Run experiment:
[tt,X,U] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, [], stop_crit, ts, maxPosErr);

%% Plot experiment:
plot_uav_exp(tt,X,U,Ts,r_margin)

function [tt,X,U] = simulate_sys(x0, xf, sim_dynamics, sim_process, controller, controller_process, supervisory_controller, stop_criterion, ts, maxPosErr,z_land)
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
            if x(3,i) >= z_land
                p_d = (xf(1:3,i)-x(1:3,i))*maxPosErr + x(1:3,i);
                x_d = [p_d;xf(4:end,i)];
                u = controller(x(:,i),x_d);
                u = controller_process(u);
                u_asif = supervisory_controller(x,u,i);
                xdot = @(t,x) sim_dynamics(x,u_asif);
                [~, x_tmp] = ode45(xdot,[t,t+ts],x(:,i));
                x(:,i) = x_tmp(end,:)';
                x(:,i) = sim_process(x(:,i), ts);
            else
                x(:,i) = taxi_uav(x(:,i),i);
            end
            X{i} = [X{i};x(:,i)'];
            U{i} = [U{i};u'];
            disp(x(1:3,:))
            %disp(u)
        end
        t = t + ts;
        tt = [tt;t];
    end
end

function plot_uav_exp(tt,X,U,Ts, r_margin)
    n_agents = size(X,2);
    [x_s, y_s, z_s] = sphere;
    for j = 1 : n_agents
       t_int = tt(1):Ts:tt(end);
       X{j} = interp1(tt(1:end-1),X{j},t_int);
    end
    n_data = length(t_int);

    %global Ts am T_exp obs r
    figure('Name','uav simulation','Position', [10 10 1800 900])
    for i=1:n_data
        clf
        hold on
        grid on
        for j = 1 : n_agents
            x = X{j}(i,:);
            plotTransforms(x(1:3), x(4:7), 'MeshFilePath', 'multirotor.stl', 'MeshColor','k','FrameSize',0.3)
            plot3(x(1),x(2),x(3),'rO');
            plot3(X{j}(1:i,1),X{j}(1:i,2),X{j}(1:i,3))
            sphere_surface = surf(x(1)+x_s*r_margin,x(2)+y_s*r_margin,x(3)+z_s*r_margin);
            if dist(X{1}(i,1:3)', X{2}(i,1:3)') <= r_margin
                 set(sphere_surface,'FaceColor',[1 0 0], ...
                'FaceAlpha',0.3,'FaceLighting','flat','EdgeAlpha',0.3)
            else
                set(sphere_surface,'FaceColor',[0 0 0], ...
                'FaceAlpha',0.1,'FaceLighting','flat','EdgeAlpha',0.3)
            end
            axis([-1 1 -1 1 1 3])
            view([0.5 2 1])
        end
        F(i) = getframe(gcf);
        drawnow
    end
    
    % Save video of experiment:
    writerObj = VideoWriter('data/uav_collision_eul.mp4','MPEG-4');
    writerObj.FrameRate = 1/Ts;
    open(writerObj);
    for i=1:length(F)
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    close(writerObj);
end
function d = dist(p1,p2)
    d = (p1-p2)'*(p1-p2);
end

function x_out = taxi_uav(x,agent_ind, h_land, ts)
    xy_park = [0.5 0 -0.5; 0 0.5 0]; 
    v = 1;
    
    v_vec = (xy_park(:,agent_ind)-x(1:2))/norm(xy_park(:,agent_ind)-x(1:2));
    xy_new = x(1:2)+v_vec*v*ts;
    x_out = [xy_new; h_land; zeros(12,1)];
end