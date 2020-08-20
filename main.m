% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

%State is defined as x = [X,Y,v,theta], u = [a,r]
clc; clear; addpath('utils')

%% Define experiment parameters:

%State constraints and backup controller parameters:
system_dynamics = @(x) dubin(x);                % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
global Ts vm rm am x_bdry
Ts = 0.1;                                       % Sampling interval
vm = 4;                                         % Maximum velocity
rm = 1.2;                                       % Maximum yaw rate
am = 2;                                         % Maximum acceleration
x_bdry = [-5,5;-5 5;vm vm;0 2*pi];              % State constraints
N_max = ceil(vm/am/Ts);                         % Maximum backup controller horizon
con1 = @(x)[-am*sign(x(3));rm];                 % Backup controller (max brake and max turn)
stop_crit1 = @(t,x)(abs(x(3))<=0);              % Stop if velocity is zero

%Koopman learning parameters:
dubin_dictionary;                               % Generate dictionary for Dubin's car system
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));  % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]
n_samples = 500;                                % Number of initial conditions to sample for training

%Collision avoidance experiment parameters:
x0 = [0;0;vm;0];
x_des = [10;0;0;0];
mpc_horizon = 20;
legacy_controller = @(x) MIQP_MPC_v3(x,x_des,20);
options = qpOASES_options('printLevel',0);
global T_exp alpha r obs
T_exp = 10;                                     % Experiment length
alpha = 1;                                      % CBF strengthening parameter
r = 1;                                          % Radius of circular obstacle
obs = [5;0];                                    % Center of circular obstacle

%% Learn approximated discrete-time Koopman operator:

X_train = collect_data(system_dynamics, con1, stop_crit1, n_samples);  %TODO: Restructure code such that it resembles more classical simulation with dynamical system + backup controller simulated forward from initial conditions
[K, C] = edmd(X_train, func_dict);
K_pows = precalc_matrix_powers(N_max,K);

%Calculate error bound (TODO)
    %-Calculate Lipschitz constant
    %-Calculate maximum residual
%Plot results (TODO)
    %-Plot data distribution
    %-Plot training prediction error and residuals

%% Evaluate Koopman approximation on test data (TODO):

%Collect test data from randomly selected initial conditions
%Evaluate prediction error
%Evaluate error bounds around predicition and compare to test data
%Plot results

%% Evaluate Koopman based CBF safety filter:

supervisory_controller = @(x,u0,N) koopman_qp_cbf(x, u0, N, system_dynamics, func_dict, K_pows, C, options);
[x_rec, u_rec, u0_rec] = run_experiment(x0, system_dynamics, legacy_controller, supervisory_controller);
plot_experiment(x_rec, u_rec, u0_rec, func_dict, K_pows, C);
