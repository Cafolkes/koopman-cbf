% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('../uav_sim_ros/codegen/','../uav_sim_ros/codegen/dynamics/','dynamics', 'controllers','koopman_learning','utils')

%% Define experiment parameters:

%State constraints and backup controller parameters:
global Ts T_max x_bdry
Ts = 0.025;                                               % Sampling interval
T_max = 1;
N_max = ceil(T_max/Ts);
ts = 1e-3;                                              % Simulator time interval
x_bdry = [-1 1; -1 1; 1 3;                              % Position limits (m)
    0 0; -pi/6 pi/6; -pi/6 pi/6;                     % Attitude limits (in euler angles ZYX) (rad)
    -1 1; -1 1; -1 1;                                   % Linear velocity limits (m/s)
    -pi/12 pi/12; -pi/12 pi/12; -pi/12 pi/12;                 % Angular velocity limits (rad/s)
    450 550; 450 550; 450 550; 450 550];                % Propeller angular velocity limits (rad/s)

% Define system and dynamics:
config = quad1_constants;
KpVxy = 0.7; %0.7
KpVz = 1; %1
KpAtt = 10; %10
KdAtt = 1; %1
KpOmegaz = 2; %2
V_max = 14.8;
V_min = 0.5;
hoverT = 0.5126*V_max; %0.52
Omega_hover = 497.61*ones(4,1);
M = [KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % Backup controller parameters

affine_dynamics = @(x) UAVDynamics(x);                  % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
backup_controller = @(x) backupU(x,M);                  % Backup controller (go to hover)
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
stop_crit1 = @(t,x)(norm(x(8:13))<=5e-2);               % Stop if velocity is zero
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config);     % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                % Processing of state data while simulating
initial_condition = @() generate_initial_state_uav();

%Koopman learning parameters:
n = 17;
func_dict = @(x) uav_D_full(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),...
         x(11),x(12),x(13),x(14),x(15),x(16),x(17));         % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]


n_samples = 100;                                         % Number of initial conditions to sample for training
gather_data = false;
tune_fit = false;

%% Learn approximated discrete-time Koopman operator:

if gather_data == true
    [T_train, X_train] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
    plot_training_data(X_train,n_samples)
    
    % Process data so it only contains the states chosen for training:
    for i = 1 : length(X_train)
        X_train{i} = X_train{i}(:,1:n);
    end
    save('data/training_data_full.mat', 'T_train','X_train');
else
    load('data/training_data_full.mat');
end

[Z, ~] = lift_data(X_train,func_dict);
n_lift = size(Z,1);
A_nom = eye(n_lift);
for i = 1 : 3
    A_nom(i+1,i+8) = Ts;
end
[Z, Z_p] = lift_data(X_train,func_dict,A_nom);

%Z = Z(:,1:100);
%Z_p = Z_p(5:end,1:100);
Z_p = Z_p(5:end,:);
if tune_fit == true
    [K, obj_vals, lambda_tuned] = edmd(Z, Z_p, 'lasso', true, [],true, 3);
    save('data/lambda_tuned_full.mat', 'lambda_tuned');
else
    load('data/lambda_tuned_full.mat');
    lambda_tuned = 1e-3*ones(n_lift,1);
    [K, obj_vals, ~] = edmd(Z, Z_p, 'lasso', true, lambda_tuned,false, 0);
end
%%
K = [zeros(4,size(Z,1)); K];
K = K + A_nom;
%K(1,1) = 1;
%for i = 1 : 3
%    K(i+1,i+1) = 1;
%    K(i+1,i+8) = Ts;
%end

%% Prepare necessary matrices and calculate error bounds:

[~,~,C] = func_dict(X_train{1}(1,:));

K_pows = precalc_matrix_powers(N_max,K);

L = 0;  
%e_max = calc_max_residual(X_train, func_dict, K, C);
tt = 0:Ts:Ts*N_max;
%error_bound = @(x) koopman_error_bound(x,X_train,L,e_max,tt,K_pows,C,func_dict);
error_bound = 0;

%% Evaluate Koopman approximation on training and test data:

if gather_data == true
    [T_test, X_test] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
    plot_training_data(X_test,n_samples)
    
    for i = 1 : length(X_test)
        X_test{i} = X_test{i}(:,1:n);
    end
    save('data/test_data_full.mat', 'T_test','X_test');
else
    load('data/test_data_full.mat');
end

% Training data fit:
fprintf('Training fit: \n')
for i = 1 : 6
    fprintf('The MSE of $x_%i$ is: %.8f \n', i+3+4, obj_vals(i+4))
end

% Test data fit:
fprintf('\nTest fit: \n')
for i = 1 : 6
    fprintf('The MSE of $x_%i$ is: %.8f \n', i+3+4, obj_vals(i+4))
end

plot_fit_uav(T_train, X_train, T_test, X_test, K_pows, C, func_dict, error_bound, N_max);

save('data/uav_learned_koopman_full.mat', 'K_pows', 'C', 'N_max');


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
        plot(tt,X{i}(:,3),'g')
        plot(tt,X{i}(:,4),':r')
        plot(tt,X{i}(:,5),':b')
        plot(tt,X{i}(:,6),':g')
        plot(tt,X{i}(:,7),':y')
        if i == 3
            title('Position and attitude data')
        end
    end

    figure(2)
    for i = 1 : n_plot
        tt = 0:Ts:(size(X{i},1)-1)*Ts;
        subplot(n_rows,n_cols,i)
        hold on
        plot(tt,X{i}(:,8),'r')
        plot(tt,X{i}(:,9),'b')
        plot(tt,X{i}(:,10),'g')
        plot(tt,X{i}(:,11),':r')
        plot(tt,X{i}(:,12),':b')
        plot(tt,X{i}(:,13),':g')
        if i==3
            title('Linear and angular velocity data')
        end
    end

    figure(3)
    for i = 1 : n_plot
        tt = 0:Ts:(size(X{i},1)-1)*Ts;
        subplot(n_rows,n_cols,i)
        hold on
        plot(tt,X{i}(:,14),'r')
        plot(tt,X{i}(:,15),'b')
        plot(tt,X{i}(:,16),'g')
        plot(tt,X{i}(:,17),'y')
        if i == 3
            title('Rotor speed data')
        end
    end
end

function plot_fit_uav(T_train, X_train, T_test, X_test, K_pows, C, func_dict, error_bound, N_max)
    global Ts x_bdry
    
    X_train_hat = predict_x(X_train, K_pows, C, func_dict);
    X_test_hat = predict_x(X_test, K_pows, C, func_dict);
    
    fig = figure(4);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    
    xyz_plot_num = [1 2; 3 4; 5 6; 7 8];
    
    for i = 1 : 2
        if i == 1
            X = X_train;
            X_hat = X_train_hat;
            T = T_train;
        else
            X = X_test;
            X_hat = X_test_hat;
            T = T_test;
        end
        
        for j = 1 : 3
            subplot(4,2,xyz_plot_num(j,i))
            hold on
            for k = 1 : length(X)
                diff = (X{k}-X_hat{k})';
                T{k} = 0 : Ts : (size(diff,2)-1)*Ts; %TODO: Fix data collection so that T is correct 
                plot(T{k},diff(j,:))
            end
            if j == 1
                ylabel('$x-\hat{x}$ (m)');
                if i == 1
                    title('Prediction error, training')
                else
                    title('Prediction error, test')
                end
            elseif j == 2
                ylabel('$y-\hat{y}$ (m)');
            elseif j == 3
                ylabel('$z-\hat{z}$ (m)');
            end
        end
        subplot(4,2,xyz_plot_num(4,i))
        hold on
        for k = 1 : length(X)
            diff = (X{k}-X_hat{k});
            T{k} = 0 : Ts : (size(diff,1)-1)*Ts; %TODO: Fix data collection so that T is correct 
            plot(T{k},vecnorm(diff(:,2:4),2,2))
        end
        xlabel('Time (sec)');
        ylabel('$||p-\hat{p}||$');
        
    end
    
    saveas(fig,'figures/uav_fit.png') 
end

function X_hat = predict_x(X, K_pows, C, func_dict)
    %Calculate fit:
    for i = 1 : length(X)
        x0 = X{i}(1,:);
        z0 = func_dict(x0);
        x_hat = [x0'];
        for j = 1 : size(X{i},1)-1
            x_hat = [x_hat C*K_pows{j}*z0];
        end
        X_hat{i} = x_hat';
    end
end
