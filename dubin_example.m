% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('utils')

%% Define experiment parameters:

%State constraints and backup controller parameters:
global Ts vm rm am x_bdry
affine_dynamics = @(x) dubin(x);                    % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                    % State is defined as x = [X,Y,v,theta], u = [a,r]
sim_dynamics = @(x,u) dubin_sim(x,u);               % System dynamics used for simulation
sim_process = @(x,ts) dubin_sim_process(x,ts);      % Processing of states while simulating
Ts = 0.1;                                           % Sampling interval
vm = 4;                                             % Maximum velocity
rm = 1.2;                                           % Maximum yaw rate
am = 2;                                             % Maximum acceleration
x_bdry = [-5,5;-5 5;vm vm;0 2*pi];                  % State constraints
N_max = ceil(vm/am/Ts);                             % Maximum backup controller horizon
con1 = @(x)[-am*sign(x(3));rm];                     % Backup controller (max brake and max turn)
stop_crit1 = @(t,x)(abs(x(3))<=0);                  % Stop if velocity is zero

%Koopman learning parameters:
dubin_dictionary;                                   % Generate dictionary for Dubin's car system
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));      % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]
n_samples = 50;                                     % Number of initial conditions to sample for training

%Collision avoidance experiment parameters:
global T_exp alpha obs r
x0 = [0;0;vm;0];                                    % Initial condition for experiment
x_des = [10;0;0;0];                                 % Desired state for legacy controller
mpc_horizon = 20;                                   % Time horizon of legacy controller
legacy_controller = @(x) MIQP_MPC_v3(x,x_des,20);   % Legacy controller (MPC)
options = qpOASES_options('printLevel',0);          % Solver options for supervisory controller
T_exp = 10;                                         % Experiment length
alpha = 1;                                          % CBF strengthening parameter
obs = [5;0];                                        % Center of circular obstacle
r = 1;                                              % Radius of circular obstacle
barrier_func = @(x) round_obs(x,obs,r);             % Barrier function

%% Learn approximated discrete-time Koopman operator:

X_train = collect_data(sim_dynamics, sim_process, con1, stop_crit1, n_samples); 
[K, C] = edmd(X_train, func_dict);
K_pows = precalc_matrix_powers(N_max,K);

% TODO: Process data so that angle always within [0,2pi]
% TODO:     Improve Lipschitz constant estimation
 
L = calc_lipschitz(4,2, affine_dynamics, con1); 
e_max = calc_max_residual(X_train, func_dict, K, C);
tt = 0:Ts:Ts*N_max;
error_bound = @(x) koopman_error_bound(x,X_train,L,e_max,tt,K_pows,C,func_dict);

plot_training_fit(X_train, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman approximation on test data:

X_test = collect_data(sim_dynamics, sim_process, con1, stop_crit1, n_samples); 
plot_test_fit(X_train, X_test, K_pows, C, func_dict, error_bound);

%% Evaluate Koopman based CBF safety filter:

supervisory_controller = @(x,u0,N) koopman_qp_cbf(x, u0, N, affine_dynamics, barrier_func, func_dict, K_pows, C, options); 
[x_rec, u_rec, u0_rec] = run_experiment(x0, sim_dynamics, sim_process, legacy_controller, supervisory_controller); 
plot_experiment(x_rec, u_rec, u0_rec, func_dict, K_pows, C);
