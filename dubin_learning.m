% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('controllers','dynamics','koopman_learning','utils','~/Documents/MATLAB/casadi-osx-matlabR2015a-v3.5.5/')

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
backup_controller = @(x) [-am*sign(x(3));rm];       % Backup controller (max brake and max turn rate)
backup_dynamics = @(x) cl_dynamics(x, affine_dynamics, backup_controller);
controller_process = @(u) u;
stop_crit1 = @(t,x)(abs(x(3))<=0);                  % Stop if velocity is zero
sim_dynamics = @(x,u) dubin_sim(x,u);               % System dynamics used for simulation
sim_process = @(x,ts) dubin_sim_process(x,ts);      % Processing of states while simulating
initial_condition = @() x_bdry(:,1)+...
    (x_bdry(:,2)-x_bdry(:,1)).*rand(4,1);           % Sample random initial value of x inside x_bdry

%Koopman learning parameters:
%dubin_dictionary;                                  % Generate dictionary for Dubin's car system
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));      % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]
n_samples = 100;                                    % Number of initial conditions to sample for training
gather_data = false;
tune_fit = false;
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
    plot_training_data(X_train,n_samples)

    save(['data/' fname '_train_data.mat'], 'T_train','X_train');
else
    load(['data/' fname '_train_data.mat']);
end

[Z, Z_p] = lift_data(X_train, func_dict);
Z_p = Z_p - Z;
Z_p = Z_p(2:end,:);
if tune_fit == true
    [K, obj_vals, lambda_tuned] = edmd(Z, Z_p, 'lasso', true, [], true, 5);
    save(['data/' fname 'lambda_tuned.mat'], 'lambda_tuned');
else
    load(['data/' fname 'lambda_tuned.mat']);
    [K, obj_vals, ~] = edmd(Z, Z_p, 'gurobi', true, lambda_tuned, false, 0);
end
K = [zeros(1,size(Z,1)); K];
K = K + eye(size(K,1));

%% Prepare necessary matrices and calculate error bounds:
[~,~,C] = func_dict(X_train{1}(1,:));
[K_pows, CK_pows] = precalc_matrix_powers(N_max,K,C);

%L = calc_lipschitz(4,2, affine_dynamics, con1); 
L = 0;  
e_max = calc_max_residual(X_train, func_dict, K, C);
tt = 0:Ts:Ts*N_max;
error_bound = @(x) koopman_error_bound(x,X_train,L,e_max,tt,K_pows,C,func_dict);

plot_training_fit(X_train, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman approximation on test data:

[T_test, X_test] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
plot_test_fit(X_train, X_test, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman based CBF safety filter:

supervisory_controller = @(x, u0, N) koopman_qp_cbf_obs(x, u0, N, affine_dynamics, backup_dynamics, barrier_func, alpha, func_dict, cell2mat(CK_pows'), options, u_lim, 4, 2);
[x_rec, u_rec, u0_rec, comp_t_rec] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller); 
plot_experiment(x_rec, u_rec, u0_rec, func_dict, CK_pows);

fprintf('\nKoopman supervisory controller avg comp. time %.2f ms, std comp. time %.2f ms\n', mean(comp_t_rec*1e3), std(comp_t_rec*1e3))

%% Evaluate integration based CBF safety filter with ODE45 (benchmark):
backup_dynamics_t = @(t,x) backup_dynamics(x);
J_cl = get_jacobian_cl(backup_dynamics, 4, 2);
sensitivity_dynamics_sim = @(t,w) sensitivity_dynamics(w, J_cl, backup_dynamics, 4);
supervisory_controller_int = @(x, u0, N) qp_cbf_obs(x, u0, N, affine_dynamics, backup_dynamics_t, barrier_func, alpha, sensitivity_dynamics_sim, options, u_lim, 4, 2);
[x_rec_int, u_rec_int, u0_rec_int, comp_t_rec_int] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller_int); 
plot_experiment_int(x_rec_int, u_rec_int, u0_rec_int, backup_dynamics_t);

fprintf('Integration based supervisory controller avg comp. time %.2f ms, std comp. time %.2f ms\n', mean(comp_t_rec_int*1e3), std(comp_t_rec_int*1e3))

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
[x_rec_int, u_rec_int, u0_rec_int, comp_t_rec_int] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller_cas); 
plot_experiment_int(x_rec_int, u_rec_int, u0_rec_int, backup_dynamics_t);

fprintf('Integration based supervisory controller avg comp. time %.2f ms, std comp. time %.2f ms\n', mean(comp_t_rec_int*1e3), std(comp_t_rec_int*1e3))

%% Save learned matrices to use in other experiments:
%close all;
save('data/dubin_learned_koopman.mat','K_pows','C','N_max')

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
