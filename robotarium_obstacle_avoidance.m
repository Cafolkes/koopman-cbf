% Initializing the agents to random positions with barrier certificates. 
% This script shows how to initialize robots to a particular point.
% Sean Wilson
% 07/2019
close all; clc; addpath('controllers','dynamics','koopman_learning','utils','data','utils/qpOASES-3.1.0/interfaces/matlab');

% Experiment parameters:
n = 3;                                                  % Number of states (unicycle model)
N = 1;                                                  % Number of robots
n_passes = 2;                                           % Number of times to travel between waypoints
radius_waypoints = 0.8;                                 % Radius of circle where waypoints are placed
start_angles = linspace(0,2*pi-(2*pi/N),N);             % Angle of ray where each robot is placed initially
end_angles = wrapTo2Pi(start_angles + pi);              % Angle of ray where each robot has final position
koopman_file = 'data/dubin_learned_koopman.mat';        % File containing learned Koopman model
r_margin = 0.08;                                        % Minimum distance between robot center points                
alpha = 1;                                              % CBF strengthening term
obs = [0;0];                                            % Center of obstacle
r_obs = 0.2;                                            % Radius of obstacle (- r_margin)

% Data storage:
file_name = 'data/obstacle_avoidance.mat';              % File to save data matrices
for i = 1 : N
    x_data{i} = [];                                     % Store state of each robot
    backup_data{i} = [];                                % Store difference between legacy and supervisory controller (norm(u0-u))
    x_init_data{i} = [];                                % Store initial point
    x_final_data{i} = [];                               % Store final point
end

% Generate intial and final positions:
initial_positions = zeros(n,N);
initial_positions(1,:) = radius_waypoints*cos(start_angles);
initial_positions(2,:) = radius_waypoints*sin(start_angles);
initial_positions(3,:) = start_angles-pi;

final_positions = zeros(n,N);
final_positions(1,:) = radius_waypoints*cos(end_angles);
final_positions(2,:) = radius_waypoints*sin(end_angles);

% Initialize robotarium environment, dynamics, and controllers:
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', initial_positions);    
si_to_uni_dynamics = create_si_to_uni_dynamics();       % Transform single integrator to unicycle dynamics
args = {'PositionError', 0.02, 'RotationError', 50};    
init_checker = create_is_initialized(args{:});          % Checks if desired position is reached
controller = create_si_position_controller();           % Legacy controller (greedy position controller)

% Construct Koopman CBF supervisory controller:
load(koopman_file)
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));
%options = qpOASES_options('printLevel',0);              % Solver options for supervisory controller
options = optimoptions('quadprog','Display','none');
affine_dynamics = @(x) dubin(x);                        % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                        % State is defined as x = [X,Y,v,theta], u = [a,r]
dt = r.time_step;
u_lim = [-r.max_linear_velocity/dt r.max_linear_velocity/dt; - pi, pi];
draw_circle(obs(1),obs(2),r_obs-r_margin);
barrier_func = @(x) round_obs(x,obs,r_obs);             % Barrier function
supervisory_controller = @(x,u0) koopman_qp_cbf_static(x, u0, N_max, affine_dynamics, barrier_func, alpha, func_dict, K_pows, C, options, u_lim);

% Get initial location data for while loop condition.
x=r.get_poses();
r.step();

v_prev = zeros(1,N);
x_prev = x;
for l = 1:n_passes
    if mod(l,2) == 0
        initial = final_positions;
        final = initial_positions;
    else
        initial = initial_positions;
        final = final_positions;
    end
    
    while(~init_checker(x, final))

        x = r.get_poses();
        dxi = controller(x(1:2, :), final(1:2, :)); 
        dxu = si_to_uni_dynamics(dxi, x);
        
        % State estimation for Dubin's car model (simplest possible
        % estimation):
        v_x_est = (x(1,:)-x_prev(1,:))/dt;
        v_y_est = (x(2,:)-x_prev(2,:))/dt;
        v_est = sqrt(v_x_est.^2+v_y_est.^2);
        x_mod = [x(1,:); x(2,:); v_est; x(3,:)]; 
        
        ad_est = (dxu(1,:) - v_est)/dt; % Estimate desired acceleration
        u_0 = [ad_est; dxu(2,:)];
        
        u_barrier = supervisory_controller(x_mod,u_0);
        vd_est = v_est + u_barrier(1)*dt;
        
        dxu = [vd_est; u_barrier(2)];
   
        r.set_velocities(1:N, dxu);
        r.step();   
        v_prev = dxu(1,:);
        x_prev = x;
        
        % Store data
        for i = 1 : N
            x_data{i} = [x_data{i} x_mod(:,i)];
            backup_data{i} = [backup_data{i} norm(u_0-u_barrier)];
            x_init_data{i} = [x_init_data{i} initial(:,i)];
            x_final_data{i} = [x_final_data{i} final(:,i)];
        end
    end    
end
save(file_name,'x_data','backup_data','x_init_data','x_final_data','obs','r_obs','r_margin');

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();