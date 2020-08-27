% Code for Koopman operator learning and CBF based safety filtering
% Written by Yuxiao Chen and Carl Folkestad
% California Institute of Technology, 2020

clc; clear; clf; close all; addpath('../uav_sim_ros/codegen/','../uav_sim_ros/codegen/dynamics/','dynamics', 'controllers','koopman_learning','utils')

%% Define experiment parameters:

%State constraints and backup controller parameters:
global Ts T_max x_bdry
Ts = 0.01;                                               % Sampling interval
T_max = 0.5;
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
KpAtt = 5; %10
KdAtt = 1; %1
KpOmegaz = 2; %2
V_max = 14.8;
V_min = 0.5;
hoverT = 0.5126*V_max; %0.52
M = [KpVxy; KpVz; KpAtt; KdAtt; KpOmegaz; hoverT];      % Backup controller parameters

affine_dynamics = @(x) UAVDynamics(x);                  % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
backup_controller = @(x) backupU(x,M);                  % Backup controller (go to hover)
controller_process = @(u) min(max(real(u),V_min*ones(4,1)),V_max*ones(4,1));
stop_crit1 = @(t,x)(norm(x(8:10))<=5e-3);               % Stop if velocity is zero
sim_dynamics = @(x,u) sim_uav_dynamics(x,u,config);     % Closed loop dynamics under backup controller
sim_process = @(x,ts) x;                                % Processing of state data while simulating
initial_condition = @() generate_initial_state_uav();

%Koopman learning parameters:
if ~exist('uav_D.m', 'file')
    uav_dictionary;                                     % Generate dictionary for uav
end
func_dict = @(x) uav_D(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),...
    x(11),x(12),x(13),x(14),x(15),x(16),x(17));         % Function dictionary, returns [D,J] = [dictionary, jacobian of dictionary]
n_samples = 100;                                         % Number of initial conditions to sample for training


%% Learn approximated discrete-time Koopman operator:

gather_data = false;
if gather_data == true
    [T_train, X_train] = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
    save('training_data.mat', 'T_train','X_train');
    plot_training_data(X_train,n_samples)
else
    load('training_data.mat');
end

[K, C] = edmd(X_train, func_dict);
K_pows = precalc_matrix_powers(N_max,K);

L = 0;  
e_max = calc_max_residual(X_train, func_dict, K, C);
tt = 0:Ts:Ts*N_max;
error_bound = @(x) koopman_error_bound(x,X_train,L,e_max,tt,K_pows,C,func_dict);
%%
plot_training_fit_uav(T_train, X_train, K_pows, C, func_dict, error_bound, N_max);

%% Evaluate Koopman approximation on test data:

%X_test = collect_data(sim_dynamics, sim_process, backup_controller, controller_process, stop_crit1, initial_condition, n_samples, ts); 
%plot_test_fit(X_train, X_test, K_pows, C, func_dict, error_bound);

function K = uav_learning(X,y)
    K_learned = zeros(n_lift);

    K_learned(1,:) = [1 zeros(1,n_lift-1)];
    for i = 2:4
        K_learned(i,:) = [zeros(1,i-1) 1 zeros(1,6) Ts zeros(1,n_lift-i-7)];
    end
    for i = 5 : n_lift
        x = Z';
        y = Z_p(i,:);
        K_learned(i,:) = lasso(x, y,'Lambda',1e-4); 
        fprintf('Learned %i out of %i predictors\n', i, n_lift);
    end
end

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

function plot_training_fit_uav(T, X, K_pows, C, func_dict, error_bound, N_max)
    global Ts x_bdry
    
    %Calculate fit:
    for i = 1 : length(X)
        x0 = X{i}(1,:);
        [z0,~] = func_dict(x0);
        x_hat = [x0'];
        for j = 1 : size(X{i},1)-1
            x_hat = [x_hat C*K_pows{j}*z0];
        end
        X_hat{i} = x_hat';
    end
        
    fig = figure(1);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    
    
%     subplot(1,2,1)
%     ind_x = 3;
%     ind_y = 4;
%     hold on
%     for i = 1 : length(X_hat)
%         X{i}(:,4) = wrapTo2Pi(X{i}(:,4));
%         scatter(X{i}(:,ind_x),X{i}(:,ind_y))
%     end
%     xlabel('Velocity ($x_3$)');
%     ylabel('Angle ($x_4$)');
%     xlim([0 x_bdry(ind_x,2)]);
%     ylim([x_bdry(ind_y,1) x_bdry(ind_y,2)]);
%     title('Training data (states)')
    
    
    err_bnd = error_bound(X{1}(1,:)');
    bound = [0];
    for k = 1 : N_max
        bound = [bound err_bnd{k}];
    end
    
    subplot(1,2,2)    
    hold on
    tt = 0 : Ts : Ts*(N_max);
    %plot(tt, bound, '--r', 'lineWidth',2)
    for i = 1 : length(X_hat)
        tt = 0 : Ts : Ts*(size(X{i},1)-1);
        diff = X{i}-X_hat{i};
        plot(tt, vecnorm(diff(:,2:4),2,2))
    end
    xlabel('Time (sec)');
    ylabel('Prediction error $||x-\hat{x}||$');
    title('Training error, position states')
    legend('Error bound')
    
    saveas(fig,'figures/uav_training_fit.png') 
    
    diff = [];
    for i = 1 : length(X_hat)
        diff = [diff X{i}'-X_hat{i}'];
    end
    diff = diff.^2;
    n_data_pts = size(diff,2);
    se = sum(diff,2)/n_data_pts;
    disp(se);
end
