
addpath('dynamics', 'controllers','koopman_learning','utils','utils/qpOASES-3.1.0/interfaces/matlab/')
file_name = 'data/uav_collision_avoidance_eul.mat';               % File to save data matrices
N = 3;
ts = 5e-3;
Ts = 1e-2;
n = 16;
m = 4;

% Define system and dynamics:
config = quad1_constants;
Kpxy = 1; %2
Kpz = 0.7; %4
KpVxy = 0.5; %0.7
KpVz = 2; %1
KpAtt = 10; %10
KdAtt = 1; %1
KpOmegaz = 2; %2
V_max = 14.8;
V_min = 0.5;
u_lim = [V_min*ones(4,1) V_max*ones(4,1)];
hoverT = 0.5126*V_max; %0.52
Omega_hover = 497.61*ones(4,1);
maxZerr = 0.1;
M = [Kpxy; Kpz; KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];           % PD controller parameters
M_backup = [KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];               % Backup controller parameters

z_land = 0.04;
stop_crit = @(x) all(x(3,:) <= z_land*ones(1,N));
legacy_controller = @(x, x_f) pdU_eul(x,x_f,M); 
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config,false,true);               % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                                % Processing of state data while simulating

r0 = 0.15;
x0_1 = [[r0*sin(2*pi/3); r0*cos(2*pi/3); 1]; zeros(9,1); Omega_hover];
x0_2 = [[(r0)*sin(0); (r0)*cos(0); 1.5]; zeros(9,1); Omega_hover];
x0_3 = [[(r0)*sin(4*pi/3); (r0)*cos(4*pi/3); 2]; zeros(9,1); Omega_hover];

xf_1 = [[0; 0; z_land]; zeros(3,1); [0;0;-2]; zeros(7,1)];
xf_2 = [[0; 0; z_land]; zeros(3,1); [0;0;-2]; zeros(7,1)];
xf_3 = [[0; 0; z_land]; zeros(3,1); [0;0;-2]; zeros(7,1)];

x0 = [x0_1 x0_2 x0_3];
xf = [xf_1 xf_2 xf_3];

% Define Koopman supervisory controller:
koopman_file = 'data/uav_learned_koopman_eul.mat';                      % File containing learned Koopman model
koopman_file_ge = 'data/uav_ge_learned_koopman.mat';                    % File containing learned Koopman model
r_margin = 0.15;                                                        % Minimum distance between robot center points                                                                              % CBF strengthening term

% Define filenames for storing data and plots:
fname_no_g_cbf = 'uav_exp_no_ground_cbf';
fname_g_cbf = 'uav_exp_ground_cbf';
fname_g_cbf_ge = 'uav_exp_ground_cbf_ge';

load(koopman_file_ge);
CK_pows_ge = CK_pows; C_ge = C; K_pows_ge = K_pows; N_max_ge = N_max;
load(koopman_file)
func_dict = @(x) uav_D_eul(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16));
func_dict_ge = @(x) uav_D_eul_ge(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16));
options = optimoptions('quadprog','Display','none');                            % Solver options for supervisory controller
affine_dynamics = @(x) UAVDynamics_eul(x);                                      % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                                                % State is defined as x = [p,q,v,w,Omega], u = [V1,V2,V3,V4]
backup_controller = @(x) backupU_eul(x,M_backup);
backup_controller_process = @(u) min(max(u,V_min*ones(4,1)),V_max*ones(4,1));
backup_dynamics = @(x) cl_dynamics(x, affine_dynamics, backup_controller, backup_controller_process); 

barrier_func_coll = @(x1,x2) collision_avoidance_3d(x1,x2,r_margin);
barrier_func_obs = @(x) paraboloid(x,5,z_land);                    % Barrier function ground

alpha = 0.25; 
supervisory_ctrl_g_cbf_ge = @(x,u0,agent_ind) koopman_qp_cbf_multi_obs_coll(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func_coll, barrier_func_obs, alpha, N, func_dict_ge, cell2mat(CK_pows_ge'), options,u_lim,16,4);

%% Define supervisory controller CBF with no ground avoidance experiment (CasADi):
import casadi.*

x = MX.sym('x', n);
q = MX.sym('q', n^2);
w = [x; q];
f_cl = sim_dynamics(x, backup_controller_process(backup_controller(x)));
J_sym = jacobian(f_cl, x);

rhs = sensitivity_dynamics_casadi(w, J_sym, f_cl, n);
ode = struct; 
ode.x = w;
ode.ode = rhs;
F = integrator('F', 'rk', ode, struct('grid', [0:Ts:N_max*Ts]));
alpha = 25;

supervisory_ctrl_no_g_cbf = @(x, u0, agent_ind) qp_cbf_multi_coll_cas(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func_coll, alpha, N, F, options, u_lim, n, m);
supervisory_ctrl_g_cbf = @(x,u0,agent_ind) qp_cbf_multi_obs_coll_cas(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func_coll, barrier_func_obs, alpha, N, F, options,u_lim,16,4);
%% Run experiments:
run_experiments = true;
if run_experiments == true
    % Run experiment with no ground CBF:
    [tt_no_g_cbf,X_no_g_cbf,U_no_g_cbf] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_ctrl_no_g_cbf, stop_crit, ts, maxZerr, z_land);

    % Run experiment with ground CBF and Koopman model trained on data with no ground effect:
    [tt_g_cbf,X_g_cbf,U_g_cbf] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_ctrl_g_cbf, stop_crit, ts, maxZerr, z_land);

    % Run experiment with ground CBF and Koopman model trained on data with ground effect:
    [tt_g_cbf_ge,X_g_cbf_ge,U_g_cbf_ge] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_ctrl_g_cbf_ge, stop_crit, ts, maxZerr, z_land);

    save('data/coll_exp.mat', 'tt_no_g_cbf','X_no_g_cbf','U_no_g_cbf','tt_g_cbf','X_g_cbf','U_g_cbf','tt_g_cbf_ge','X_g_cbf_ge','U_g_cbf_ge');
else
    load('data/coll_exp.mat');
end

%% Plot experiment:

plot_uav_exp(tt_no_g_cbf,X_no_g_cbf,U_no_g_cbf,1/30,r_margin,z_land, false,fname_no_g_cbf)
plot_uav_exp(tt_g_cbf,X_g_cbf,U_g_cbf,1/30,r_margin,z_land,false,fname_g_cbf)
plot_uav_exp(tt_g_cbf_ge,X_g_cbf_ge,U_g_cbf_ge,1/30,r_margin,z_land,false,fname_g_cbf_ge)

%%
plot_exp_summary(tt_no_g_cbf, X_no_g_cbf, U_no_g_cbf, tt_g_cbf, X_g_cbf, U_g_cbf, tt_g_cbf_ge, X_g_cbf_ge, U_g_cbf_ge,z_land, r_margin)

%%% Supporting functions:
function [tt,X,U] = simulate_sys(x0, xf, sim_dynamics, sim_process, controller, controller_process, supervisory_controller, stop_criterion, ts, maxZErr,z_land)
    t = 0;
    tt = 0;
    n_agents = size(x0,2);
    x = x0;
    x_prev = x0;
    landed = false(1,n_agents);
    parked = false(1,n_agents);
    for i = 1 : n_agents
        X{i} = [];
        U{i} = [];
    end
    
    while ~all(parked) && t <= 6
        for i = 1 : n_agents
            if landed(i) == false
                z_d = (xf(3,i)-x(3,i))*maxZErr + x(3,i);
                x_d = [xf(1:2,i); z_d; xf(4:end,i)];
                u = controller(x(:,i),x_d);
                u = controller_process(u);
                u_asif = supervisory_controller(x_prev,u,i);
                xdot = @(t,x) sim_dynamics(x,u_asif);
                [~, x_tmp] = ode45(xdot,[t,t+ts],x(:,i));
                x(:,i) = x_tmp(end,:)';
                x(:,i) = sim_process(x(:,i), ts);
                
                if x(3,i) <= z_land+0.01
                    landed(i) = true;
                end
            else
                [x(:,i), parked(i)] = taxi_uav(x(:,i),i, z_land, ts);
            end
            X{i} = [X{i};x(:,i)'];
            U{i} = [U{i};u'];
            disp(x(1:3,:))
            %disp(t)
        end
        x_prev = x;
        t = t + ts;
        tt = [tt;t];
    end
end

function plot_uav_exp(tt,X,U,Ts, r_margin, z_land, draw, fname)
    n_agents = size(X,2);
    [x_s, y_s, z_s] = sphere;
    for j = 1 : n_agents
       t_int = tt(1):Ts:tt(end);
       X{j} = interp1(tt(1:end-1),X{j},t_int);
    end
    n_data = length(t_int)-1;

    %global Ts am T_exp obs r
    figure('Name','uav simulation','Position', [10 10 1200 1000])
    for i=1:n_data
        clf
        hold on
        grid on
        for j = 1 : n_agents
            x = X{j}(i,:);
            % Plot quadrotor model:
            plotTransforms(x(1:3), eul2quat(x(4:6),'XYZ'), 'MeshFilePath', 'multirotor.stl', 'MeshColor',[0.2 0.2 0.2],'FrameSize',0.3)
            
            % Plot trajectory trace:
            start_ind = max(1,i-80);
            plot3(X{j}(start_ind:i,1),X{j}(start_ind:i,2),X{j}(start_ind:i,3))
            
            % Plot sphere:
            sphere_surface(j) = surf(x(1)+x_s*r_margin,x(2)+y_s*r_margin,x(3)+z_s*r_margin);
            if dist(X{1}(i,1:3)', X{2}(i,1:3)') <= (r_margin)^2
                 set(sphere_surface(j),'FaceColor',[1 0 0], ...
                'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor',[0.1 0.1 0.1])
            elseif x(3) <= z_land+0.0025
                set(sphere_surface(j),'FaceColor',[0 1 0], ...
                'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
            else
                set(sphere_surface(j),'FaceColor',[0 0 0], ...
                'FaceAlpha',0.1,'FaceLighting','gouraud','EdgeColor',[0.8 0.8 0.8], 'EdgeAlpha',0.2)
            end
            
            % Plot landing pad:
            draw_circle_3d(0,0,0.25,'w',5)
            
            % Set plot parameters:
            axis([-0.6 0.6 -0.8 0.6 0 2.])
            view([0.5 2 1])
            box on
            
            % Set face color of xy-plane:
            XL = get(gca, 'XLim');
            YL = get(gca, 'YLim');
            patch([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], [0 0 0 0], 'FaceColor', [0.7 0.7 0.7]);
        end
        F(i) = getframe(gcf);
        if draw == true
            drawnow
        end
    end
    
    % Save video of experiment:
    writerObj = VideoWriter(['figures/' fname '.mp4'],'MPEG-4');
    writerObj.FrameRate = floor(1/Ts);
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

function plot_exp_summary(tt_no_g_cbf, X_no_g_cbf, U_no_g_cbf, tt_g_cbf, X_g_cbf, U_g_cbf, tt_g_cbf_ge, X_g_cbf_ge, U_g_cbf_ge, z_land, r_margin)
    n_agents = size(X_no_g_cbf,2);
    % Plot altitude and velocities of the 3 landing scenarios:
    figure('Name','uav experiment summary','Position',[560 528 1680 450])
    fs = 16;
    color_map = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
    lw = 2;
    z_ind = 3;
    v_z_ind = 9;
    
    for i = 1 : n_agents
        subplot(3,3,2)
        hold on
        plot(tt_no_g_cbf(2:end), X_no_g_cbf{i}(:,z_ind)-z_land,'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_no_g_cbf(2:marker_spacing:end), X_no_g_cbf{i}(1:marker_spacing:end,z_ind)-z_land, strcat('',line_specifiers(i)));
        %p(i) = plot(tt_no_g_cbf(end),X_no_g_cbf{i}(end,z_ind),strcat('-',line_specifiers(i)));
        
        subplot(3,3,3)
        hold on
        plot(tt_no_g_cbf(2:end), X_no_g_cbf{i}(:,v_z_ind),'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_no_g_cbf(2:marker_spacing:end), X_no_g_cbf{i}(1:marker_spacing:end,v_z_ind), strcat('',line_specifiers(i)));
        %p(3+i) = plot(tt_no_g_cbf(end),X_no_g_cbf{i}(end,v_z_ind),strcat(':',line_specifiers(i)));

        subplot(3,3,5)
        hold on
        plot(tt_g_cbf(2:end), X_g_cbf{i}(:,z_ind)-z_land,'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_g_cbf(2:marker_spacing:end), X_g_cbf{i}(1:marker_spacing:end,z_ind)-z_land, strcat('',line_specifiers(i)));
        
        subplot(3,3,6)
        hold on
        plot(tt_g_cbf(2:end), X_g_cbf{i}(:,v_z_ind),'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_g_cbf(2:marker_spacing:end), X_g_cbf{i}(1:marker_spacing:end,v_z_ind), strcat('',line_specifiers(i)));

        subplot(3,3,8)
        hold on
        plot(tt_g_cbf_ge(2:end), X_g_cbf_ge{i}(:,z_ind)-z_land,'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_g_cbf_ge(2:marker_spacing:end), X_g_cbf_ge{i}(1:marker_spacing:end,z_ind)-z_land, strcat('',line_specifiers(i)));
        %ylabel('Altitude (m)')
        subplot(3,3,9)
        hold on
        plot(tt_g_cbf_ge(2:end), X_g_cbf_ge{i}(:,v_z_ind),'LineWidth', lw, 'Color',color_map(i,:))
        %hold on
        %plot(tt_g_cbf_ge(2:marker_spacing:end), X_g_cbf_ge{i}(1:marker_spacing:end,v_z_ind), strcat('',line_specifiers(i)));
        
        
    end
    subplot(3,3,2)
    xlim([0 5])
    ylabel('Alt (m)')
    set(gca,'FontSize',fs)
    title(' ')
    
    subplot(3,3,3)
    xlim([0 5])
    ylabel('Vel (m/s)')
    %legend(p,'Altitude UAV 1', 'Altitude UAV 2', 'Altitude UAV 3', 'Velocity UAV 1', 'Velocity UAV 2', 'Velocity UAV 3','NumColumns',2)
    set(gca,'FontSize',fs)
    title('Scenario 1: no ground CBF, supervisory controller based on dynamics model')
    
    subplot(3,3,5)
    xlim([0 5])
    ylabel('Alt (m)')
    set(gca,'FontSize',fs)
    title(' ')
    
    subplot(3,3,6)
    xlim([0 5])
    ylabel('Vel (m/s)')
    title('Scenario 2: with ground CBF, supervisory controller based on dynamics model')
    set(gca,'FontSize',fs)
    
    subplot(3,3,8)
    xlim([0 5])
    ylabel('Alt (m)')
    xlabel('Time (sec)')
    set(gca,'FontSize',fs)
    title(' ')
    
    subplot(3,3,9)
    xlim([0 5])
    ylabel('Vel (m/s)')
    xlabel('Time (sec)')
    title('Scenario 3: with ground CBF, supervisory controller based on Koopman operator')
    set(gca,'FontSize',fs)
    
    % Plot 3D trajectories and uav snapshot (similar to video sim):
    [x_s, y_s, z_s] = sphere;
    %figure('Name','3D trajectories')
    subplot(1,3,1)
    hold on
    grid on
    
    X = X_g_cbf_ge; 
    T = tt_g_cbf_ge; 
    uav_pos_ind = 230;
    
    for i = 1 : n_agents
        x = X{i}; 
        % Plot quadrotor model:
        plotTransforms(x(uav_pos_ind,1:3), eul2quat(x(uav_pos_ind,4:6),'XYZ'), 'MeshFilePath', 'multirotor.stl', 'MeshColor',[0.2 0.2 0.2],'FrameSize',0.3)

        % Plot trajectory trace:
        p3(i) = plot3(x(:,1),x(:,2),x(:,3),'LineWidth', lw, 'Color',color_map(i,:));

        % Plot sphere:
        sphere_surface(i) = surf(x(uav_pos_ind,1)+x_s*r_margin,x(uav_pos_ind,2)+y_s*r_margin,x(uav_pos_ind,3)+z_s*r_margin);
        if dist(X{1}(uav_pos_ind,1:3)', X{2}(uav_pos_ind,1:3)') <= (r_margin)^2
             set(sphere_surface(i),'FaceColor',[1 0 0], ...
            'FaceAlpha',0.3,'FaceLighting','gouraud','EdgeColor',[0.1 0.1 0.1])
        elseif x(uav_pos_ind,3) <= z_land+0.001
            set(sphere_surface(i),'FaceColor',[0 1 0], ...
            'FaceAlpha',0.2,'FaceLighting','gouraud','EdgeColor','none')
        else
            set(sphere_surface(i),'FaceColor',[0 0 0], ...
            'FaceAlpha',0.1,'FaceLighting','gouraud','EdgeColor',[0.8 0.8 0.8], 'EdgeAlpha',0.2)
        end

        % Plot landing pad:
        draw_circle_3d(0,0,0.25,'w',5)

        % Set plot parameters:
        axis([-0.4 0.4 -0.4 0.4 0 1.5])
        view([1 2 1])
        box on

        % Set face color of xy-plane:
        XL = get(gca, 'XLim');
        YL = get(gca, 'YLim');
        patch([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], [0 0 0 0], 'FaceColor', [0.7 0.7 0.7]);
    end
    legend(p3, 'UAV 1','UAV 2','UAV 3','Location', 'northwest') 
    title('Snapshot of experiment ')
    set(gca,'FontSize',fs)
    saveas(gcf, 'figures/uav_coll_summary_raw.png')
end

function [x_out, parked] = taxi_uav(x,agent_ind, z_land, ts)
    xy_park = [0.4 0.4 0; -0.4 -0 -0.4]; 
    v = min(0.5,4*max(abs(0.5-norm(xy_park(:,agent_ind)-x(1:2))),0.1))*1;
    parked=false;
    if norm(xy_park(:,agent_ind)-x(1:2)) >= 1e-2
        v_vec = (xy_park(:,agent_ind)-x(1:2))/norm(xy_park(:,agent_ind)-x(1:2));
        xy_new = x(1:2)+v_vec*v*ts;
        x_out = [xy_new; z_land; zeros(13,1)];
    else
        x_out = x;
        parked=true;
    end
end

function [] =draw_circle_3d(X,Y,R,color,linewidth)
    if nargin<4
        color='r';

    end
    if nargin<5
        linewidth=1;
    end
    theta=(1:360)*pi/180;
    for i=1:size(theta,2)
        x(i)=X+R*cos(theta(i));
        y(i)=Y+R*sin(theta(i));
        z(i)=0.005;
    end
    plot3(x,y,z,'r','linewidth',linewidth);
    fill3(x,y,z,color);
end