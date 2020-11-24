% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('controllers','dynamics','koopman_learning','utils','~/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.5/', 'utils/qpOASES-3.1.0/interfaces/matlab/', '~/Documents/MATLAB/Intlab_V12/')

%% Define experiment parameters:

%State constraints and backup controller parameters:
global Ts T_max vm rm am x_bdry
Ts = 0.1;                                           % Sampling interval
T_max = 2;
vm = 0.15;                                          % Maximum velocity
rm = pi;                                            % Maximum yaw rate
am = 0.1;                                           % Maximum acceleration
u_lim = [-am am; -rm, rm];
x_bdry = [-1.6 1.6;-1 1;vm vm;0 2*pi];              % State constraints
ts = 0.02;
N_max = ceil(vm/am/Ts)+1;                           % Maximum backup controller horizon
n = 4;
m = 2;

affine_dynamics = @(x) dubin(x);                    % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                    % State is defined as x = [X,Y,v,theta], u = [a,r]
backup_controller = @(x) [-am*sign(x(3)); abs(x(3))*rm/vm];       % Backup controller (max brake and max turn rate)
backup_controller_process = @(u) u;
backup_dynamics = @(x) cl_dynamics(x, affine_dynamics, backup_controller, backup_controller_process);
controller_process = @(u) u;
stop_crit1 = @(t,x)(abs(x(3))<=0);                  % Stop if velocity is zero
sim_dynamics = @(x,u) dubin_sim(x,u);               % System dynamics used for simulation
sim_process = @(x,ts) dubin_sim_process(x,ts);      % Processing of states while simulating
initial_condition = @() x_bdry(:,1)+ ...
    (x_bdry(:,2)-x_bdry(:,1)).*rand(4,1);           % Sample random initial value of x inside x_bdry

%Koopman learning parameters:
dubin_dictionary(true);                             % Generate dictionary for Dubin's car system
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));      % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]
n_samples = 200;                                    % Number of initial conditions to sample for training
gather_data = false;
fname = 'dubin';

%Collision avoidance experiment parameters:
global T_exp obs r
x0 = [0;0;vm;0];                                    % Initial condition for experiment
x_des = [1;0;0;0];                                  % Desired state for legacy controller
mpc_horizon = 20;                                   % Time horizon of legacy controller
legacy_controller = @(x) MIQP_MPC_v3(x,x_des,20);   % Legacy controller (MPC)
options = optimoptions('quadprog','Display','none');% Solver options for supervisory controller
T_exp = 10;                                         % Experiment length
alpha = 2;                                          % CBF strengthening parameter
obs = [0.5;0];                                      % Center of circular obstacle
r = 0.05;                                           % Radius of circular obstacle
barrier_func = @(x) round_obs(x,obs,r);             % Barrier function

%% Learn approximated discrete-time Koopman operator:
if gather_data == true
    [T_train, X_train] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
    plot_training_data(X_train, n_samples)

    save(['data/' fname '_train_data.mat'], 'T_train','X_train');
else
    load(['data/' fname '_train_data.mat']);
end

[Z, Z_p] = lift_data(X_train, func_dict, false, true);
Z_p = Z_p - Z;

[K, obj_vals, ~] = edmd(Z, Z_p, 'gurobi', true, 1e-3, 1e-3, true, false, 0);
[~,~,C] = func_dict(X_train{1}(1,:));

K(1,:) = zeros(1,size(K,2)); % Override constant term learning
K = K + eye(size(K,1));

%% Train Koopman operator for error bound estimation:

[Z_bound, Z_p_bound] = lift_data(X_train, func_dict, true, true);
Z_p_bound = Z_p_bound-Z_bound;

[K_bound, obj_vals, ~] = edmd(Z_bound, Z_p_bound, 'gurobi', true, 1e-3, 1e-3, false, false, 0);
K_bound(1,:) = zeros(1,size(Z_bound,1)); % Override constant term learning
K_bound = K_bound + eye(size(K_bound,1));

eig_max = abs(max(eig(K_bound)));
mu_min = 0;
    
%% Prepare necessary matrices and calculate error bounds:
[K_pows, CK_pows] = precalc_matrix_powers(N_max,K,C);
[K_pows_bound, ~] = precalc_matrix_powers(N_max,K_bound,C);

func_dict = @(x) func_dict(x);
n_lift = length(func_dict(ones(4,1)));
L = calc_lipschitz(n_lift, func_dict);

e_max = calc_max_residual(X_train, func_dict, K, C);
tt = 0:Ts:Ts*N_max;
error_bound = @(x) koopman_error_bound_mu(mu_min, L,e_max,tt,K_pows,C(1:2,:),func_dict,1);

plot_training_fit(X_train, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman approximation on test data:

[T_test, X_test] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
plot_test_fit(X_train, X_test, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman based CBF safety filter:

supervisory_controller = @(x, u0, N) koopman_qp_cbf_obs(x, u0, N, affine_dynamics, backup_dynamics, barrier_func, alpha, func_dict, cell2mat(CK_pows'), options, u_lim, 4, 2);
[x_rec, u_rec, u0_rec, comp_t_rec, int_t_rec] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller); 
plot_experiment(x_rec, u_rec, u0_rec, func_dict, CK_pows);

fprintf('\nKoopman CBF supervisory controller:\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec*1e3), std(comp_t_rec*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec*1e3), std(int_t_rec*1e3))

%% Evaluate robust Koopman based CBF safety filter:

supervisory_controller_rob = @(x, u0, N) koopman_qp_cbf_obs_rob(x, u0, N, affine_dynamics, backup_dynamics, barrier_func, alpha, func_dict, cell2mat(CK_pows'), options, u_lim, 4, 2, error_bound);
[x_rec_rob, u_rec_rob, u0_rec_rob, comp_t_rec_rob, int_t_rec_rob] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller_rob); 
plot_experiment(x_rec_rob, u_rec_rob, u0_rec_rob, func_dict, CK_pows);

fprintf('\nRobust Koopman CBF supervisory controller:\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec_rob*1e3), std(comp_t_rec_rob*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec_rob*1e3), std(int_t_rec_rob*1e3))

%% Evaluate integration based CBF safety filter with ODE45 (benchmark):
backup_dynamics_t = @(t,x) backup_dynamics(x);
J_cl = get_jacobian_cl(backup_dynamics, 4, 2);
sensitivity_dynamics_sim = @(t,w) sensitivity_dynamics(w, J_cl, backup_dynamics, 4);
supervisory_controller_int = @(x, u0, N) qp_cbf_obs(x, u0, N, affine_dynamics, backup_dynamics_t, barrier_func, alpha, sensitivity_dynamics_sim, options, u_lim, 4, 2);
[x_rec_ode45, u_rec_ode45, u0_rec_ode45, comp_t_rec_ode45, int_t_rec_ode45] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller_int); 
plot_experiment_int(x_rec_ode45, u_rec_ode45, u0_rec_ode45, backup_dynamics_t);

fprintf('\nIntegration based CBF supervisory controller (ODE45):\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec_ode45*1e3), std(comp_t_rec_ode45*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec_ode45*1e3), std(int_t_rec_ode45*1e3))

%% Evaluate integration based CBF safety filter with casadi (benchmark):
import casadi.*
backup_dynamics_t = @(t,x) backup_dynamics(x);

x = MX.sym('x', n);
q = MX.sym('q', n^2);
w = [x; q];
f_cl = backup_dynamics(x);
J_sym = jacobian(f_cl, x);

rhs = sensitivity_dynamics_casadi(w, J_sym, f_cl, n);
ode = struct;
ode.x = w;
ode.ode = rhs;
F = integrator('F', 'rk', ode, struct('grid', [0:Ts:N_max*Ts]));

supervisory_controller_cas = @(x, u0, N) qp_cbf_obs_cas(x, u0, N, affine_dynamics, backup_dynamics_t, barrier_func, alpha, F, options, u_lim, 4, 2);
[x_rec_cas, u_rec_cas, u0_rec_cas, comp_t_rec_cas, int_t_rec_cas] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller_cas); 
plot_experiment_int(x_rec_cas, u_rec_cas, u0_rec_cas, backup_dynamics_t);

fprintf('\nIntegration based CBF supervisory controller (casADi):\n')
fprintf('Average computation time %.2f ms, std computation time %.2f ms\n', mean(comp_t_rec_cas*1e3), std(comp_t_rec_cas*1e3))
fprintf('Average integration time %.2f ms, std computation time %.2f ms\n', mean(int_t_rec_cas*1e3), std(int_t_rec_cas*1e3))

%% Save learned matrices to use in other experiments:

save('data/dubin_learned_koopman.mat','CK_pows','C','N_max')

%% Supporting functions:
function plot_training_data(X,n_samples)
    global Ts
    n_plot = min(n_samples,100);
    n_rows = 10;
    n_cols = 10;
    figure(1)
    for i = 1 : n_plot
        tt = 0:Ts:(size(X{i},1)-1)*Ts;
        subplot(n_rows,n_cols,i)
        hold on
        plot(tt,X{i}(:,1),'r')
        plot(tt,X{i}(:,2),'b')   
        plot(tt,X{i}(:,3),':r')
        plot(tt,X{i}(:,4),':b')
        if i == 3
            title('State data')
        end
    end
end
