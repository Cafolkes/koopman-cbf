% % 
clear; close all;
addpath('dynamics', 'controllers','koopman_learning','utils','utils/qpOASES-3.1.0/interfaces/matlab/', '~/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.5/', 'utils/qpOASES-3.1.0/interfaces/matlab/')
file_name = 'data/uav_collision_avoidance_eul.mat';               % File to save data matrices
N = 2;
n = 16;
m = 4;
ts = 1e-2;
global Ts 
Ts = 1e-2;

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
V_min = 0; % 1
u_lim = [V_min*ones(4,1) V_max*ones(4,1)];
hoverT = 0.5126*V_max; %0.52
Omega_hover = 497.61*ones(4,1);
maxPosErr = 0.3;
M = [Kpxy; Kpz; KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % PD controller parameters
M_backup = [KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % Backup controller parameters

stop_crit = @(x,x_f) norm(x(1:3,1)-x_f(1:3,1)) <= 1e-1 && norm(x(1:3,2)-x_f(1:3,2)) <= 1e-1;
legacy_controller = @(x, x_f) pdU_eul(x,x_f,M); 
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config,false,false);     
sim_process = @(x,ts) x;                                % Processing of state data while simulating

x0_1 = [[-0.5; 0.5; 1.5]; zeros(9,1); Omega_hover];
x0_2 = [[0.5; -0.5; 2]; zeros(9,1); Omega_hover];
xf_1 = [[0.5; -0.5; 2]; zeros(9,1); Omega_hover];
xf_2 = [[-0.5; 0.5; 1.5]; zeros(9,1); Omega_hover];
x0 = [x0_1 x0_2];
xf = [xf_1 xf_2];

% Define Koopman supervisory controller:
koopman_file = 'data/uav_learned_koopman_eul.mat';          % File containing learned Koopman model
r_margin = 0.15;                                            % Minimum distance between robot center points                
alpha =7.5;                                                % CBF strengthening term

load(koopman_file)
func_dict = @(x) uav_D_eul(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16));
options = optimoptions('quadprog','Display','none');                % Solver options for supervisory controller
affine_dynamics = @(x) UAVDynamics_eul(x);                          % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                                    % State is defined as x = [p,q,v,w,Omega], u = [V1,V2,V3,V4]
barrier_func = @(x1,x2) collision_avoidance_3d(x1,x2,r_margin);     % Barrier function
backup_controller = @(x) backupU_eul(x,M_backup);
backup_controller_process = @(u) min(max(u,V_min*ones(4,1)),V_max*ones(4,1));
backup_dynamics = @(x) cl_dynamics(x, affine_dynamics, backup_controller, backup_controller_process); 
supervisory_controller = @(x,u0,agent_ind) koopman_qp_cbf_multi_coll(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func, alpha, N, func_dict, cell2mat(CK_pows'), options, u_lim,16,4);

%% Run experiment with Koopman CBF safety filter:
[tt,X,U, comp_t_rec, int_t_rec] = simulate_sys(x0,xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_controller, stop_crit, ts, maxPosErr);

fprintf('\nKoopman CBF supervisory controller:\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec*1e3), std(comp_t_rec*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec*1e3), std(int_t_rec*1e3))

%% Evaluate integration based CBF safety filter with ODE45 (benchmark):
x_sym = sym('x_sym',[16,1],'real');
f_cl = sim_dynamics(x_sym, backup_controller(x_sym));
J_sym = jacobian(f_cl, x_sym);
J_cl = matlabFunction(J_sym, 'Vars', {x_sym});

f_cl_sim = @(x) sim_dynamics(x, backup_controller_process(backup_controller(x)));
sensitivity_dynamics_sim = @(t,w) sensitivity_dynamics(w, J_cl, f_cl_sim, n);
supervisory_controller_ode45 = @(x, u0, agent_ind) qp_cbf_multi_coll(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func, alpha, N, sensitivity_dynamics_sim, options, u_lim, n, m);

[tt_ode45, X_ode45, U_ode45, comp_t_rec_ode45, int_t_rec_ode45] = simulate_sys(x0, xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_controller_ode45, stop_crit, ts, maxPosErr);

fprintf('\nIntegration based CBF supervisory controller (ODE45):\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec_ode45*1e3), std(comp_t_rec_ode45*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec_ode45*1e3), std(int_t_rec_ode45*1e3))

%% Evaluate integration based CBF safety filter with casadi (benchmark):
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

supervisory_controller_cas = @(x, u0, agent_ind) qp_cbf_multi_coll_cas(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func, alpha, N, F, options, u_lim, n, m);
[tt_cas, X_cas, U_cas, comp_t_rec_cas, int_t_rec_cas] = simulate_sys(x0, xf, sim_dynamics, sim_process, legacy_controller, controller_process, supervisory_controller_cas, stop_crit, ts, maxPosErr);

fprintf('\nIntegration based CBF supervisory controller (casADi):\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec_cas*1e3), std(comp_t_rec_cas*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec_cas*1e3), std(int_t_rec_cas*1e3))

%% Plot experiment:
plot_uav_exp(tt,X,U,1/24,r_margin, 'koop')  % Plot Koopman CBF experiment
plot_uav_exp(tt_ode45,X_ode45,U_ode45,1/24,r_margin, 'ode45')  % Plot ode45 CBF experiment
plot_uav_exp(tt_cas,X_cas,U_cas,1/24,r_margin, 'casadi')  % Plot casADi CBF experiment


function [tt,X,U, comp_t_rec, int_t_rec] = simulate_sys(x0, xf, sim_dynamics, sim_process, controller, controller_process, supervisory_controller, stop_criterion, ts, maxPosErr)
    t = 0;
    tt = 0;
    n_agents = size(x0,2);
    x = x0;
    comp_t_rec = [];
    int_t_rec = [];
    for i = 1 : n_agents
        X{i} = [];
        U{i} = [];
    end
    
    while ~stop_criterion(x,xf)
        for i = 1 : n_agents
            p_d = (xf(1:3,i)-x(1:3,i))*maxPosErr + x(1:3,i);
            x_d = [p_d;xf(4:end,i)];
            x_cur = x;
            u = controller(x_cur(:,i),x_d);
            u = controller_process(u);
            comp_t0 = posixtime(datetime('now'));
            [u_asif, int_time] = supervisory_controller(x_cur,u,i);
            comp_tf = posixtime(datetime('now')) - comp_t0;
            xdot = @(t,x) sim_dynamics(x,u_asif);
            [~, x_tmp] = ode45(xdot,[t,t+ts],x_cur(:,i));
            x(:,i) = x_tmp(end,:)';
            x(:,i) = sim_process(x(:,i), ts);
            X{i} = [X{i};x(:,i)'];
            U{i} = [U{i};u'];
            comp_t_rec = [comp_t_rec comp_tf];
            int_t_rec = [int_t_rec int_time];
        end
        t = t + ts;
        tt = [tt;t];
    end
end

function plot_uav_exp(tt,X,U,Ts, r_margin, fname)
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
            plotTransforms(x(1:3), eul2quat(x(4:6),'XYZ'), 'MeshFilePath', 'multirotor.stl', 'MeshColor',[0.2 0.2 0.2],'FrameSize',0.3)
            plot3(x(1),x(2),x(3),'rO');
            start_ind = max(1,i-80);
            plot3(X{j}(start_ind:i,1),X{j}(start_ind:i,2),X{j}(start_ind:i,3))
            sphere_surface = surf(x(1)+x_s*r_margin,x(2)+y_s*r_margin,x(3)+z_s*r_margin);
            if dist(X{1}(i,1:3)', X{2}(i,1:3)') <= (2*r_margin)^2
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
    writerObj = VideoWriter(['figures/uav_collision_' fname '.mp4'],'MPEG-4');
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