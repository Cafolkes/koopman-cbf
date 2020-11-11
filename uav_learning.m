% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('../uav_sim_ros/codegen/','../uav_sim_ros/codegen/dynamics/','dynamics', 'controllers','koopman_learning','utils')

%% Define experiment parameters:

%State constraints and backup controller parameters:
global Ts T_max x_bdry
Ts = 0.01;                                               % Sampling interval
T_max = 1;
N_max = ceil(T_max/Ts);
ts = 1e-3;                                              % Simulator time interval
x_bdry = [-1 1; -1 1; 0.2 2;                              % Position limits (m)
    -pi/6 pi/6; -pi/6 pi/6; -pi/12 pi/12;                     % Attitude limits (in euler angles XYZ) (rad)
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
z_land = 0.05;

affine_dynamics = @(x) UAVDynamics_eul(x);                  % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
backup_controller = @(x) backupU_eul(x,M);                  % Backup controller (go to hover)
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
stop_crit1 = @(t,x)(norm(x(7:12))<=5e-2 || x(3) <= z_land);               % Stop if velocity is zero
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config,false,false);     % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                % Processing of state data while simulating
initial_condition = @() generate_initial_state_uav(false);
fname = 'uav';

%Koopman learning parameters:
n = 16;
func_dict = @(x) uav_D_eul_ge(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),...
         x(11),x(12),x(13),x(14),x(15),x(16));         % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]

n_samples = 250;                                         % Number of initial conditions to sample for training
gather_data = true;
tune_fit = true;

%% Learn approximated discrete-time Koopman operator:

if gather_data == true
    [T_train, X_train] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
    plot_training_data(X_train,n_samples)
    
    % Process data so it only contains the states chosen for training:
    for i = 1 : length(X_train)
        X_train{i} = X_train{i}(:,1:n);
    end
    save(['data/' fname '_train_data.mat'], 'T_train','X_train');
else
    load(['data/' fname '_train_data.mat']);
end

[Z, Z_p] = lift_data(X_train,func_dict);
Z_p = Z_p - Z;
Z_p = Z_p(5:end,:);
if tune_fit == true
    [K, obj_vals, lambda_tuned] = edmd(Z, Z_p, 'lasso', true, [],true, 5);
    save(['data/' fname '_lambda_tuned.mat'], 'lambda_tuned');
else
    load(['data/' fname '_lambda_tuned.mat']);
    [K, obj_vals, ~] = edmd(Z, Z_p, 'lasso', true, lambda_tuned,false, 0);
    %[K, obj_vals, ~] = edmd(Z, Z_p, 'gurobi', true, lambda_tuned, false, 0);
end
%%
K = [zeros(4,size(Z,1)); K];
K = K + eye(size(K,1));
for i = 1 : 3
    K(i+1,i+7) = Ts;
end

%% Prepare necessary matrices and calculate error bounds:

[~,~,C] = func_dict(X_train{1}(1,:));

[K_pows, CK_pows] = precalc_matrix_powers(N_max,K,C);

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
    save(['data/' fname '_test_data.mat'], 'T_test','X_test');
else
    load(['data/' fname '_test_data.mat']);
end

% Training data fit:
fprintf('Training fit: \n')
for i = 1 : 6
    fprintf('The MSE of $x_%i$ is: %.8f \n', i+6, obj_vals(i+3))
end

% Test data fit:
fprintf('\nTest fit: \n')
for i = 1 : 6
    fprintf('The MSE of $x_%i$ is: %.8f \n', i+6, obj_vals(i+3))
end

plot_fit_uav(T_train, X_train, T_test, X_test, K_pows, C, func_dict, error_bound, N_max, fname);

save(['data/' fname '_learned_koopman.mat'], 'K_pows', 'CK_pows', 'C', 'N_max');


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
        if i == 3
            title('Position and attitude data')
        end
    end

    figure(2)
    for i = 1 : n_plot
        tt = 0:Ts:(size(X{i},1)-1)*Ts;
        subplot(n_rows,n_cols,i)
        hold on
        plot(tt,X{i}(:,7),'r')
        plot(tt,X{i}(:,8),'b')
        plot(tt,X{i}(:,9),'g')
        plot(tt,X{i}(:,10),':r')
        plot(tt,X{i}(:,11),':b')
        plot(tt,X{i}(:,12),':g')
        if i==3
            title('Linear and angular velocity data')
        end
    end

    figure(3)
    for i = 1 : n_plot
        tt = 0:Ts:(size(X{i},1)-1)*Ts;
        subplot(n_rows,n_cols,i)
        hold on
        plot(tt,X{i}(:,13),'r')
        plot(tt,X{i}(:,14),'b')
        plot(tt,X{i}(:,15),'g')
        plot(tt,X{i}(:,16),'y')
        if i == 3
            title('Rotor speed data')
        end
    end
end

function plot_fit_uav(T_train, X_train, T_test, X_test, K_pows, C, func_dict, error_bound, N_max, fname)
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
    
    saveas(fig,['figures/' fname '_fit.png']) 
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
