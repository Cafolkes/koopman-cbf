% Initializing the agents to random positions with barrier certificates. 
% This script shows how to initialize robots to a particular point.
% Sean Wilson
% 07/2019
close all; clc; addpath('controllers','dynamics','koopman_learning','utils');

% Experiment parameters:
n = 3;                                              % Number of states (unicycle model)
N = 5;                                              % Number of robots
n_passes = 2;                                       % Number of times to travel between waypoints
radius_waypoints = 0.8;                             % Radius of circle where waypoints are placed
start_angles = linspace(0,2*pi-(2*pi/N),N);         % Angle of ray where each robot is placed initially
end_angles = wrapTo2Pi(start_angles + pi);          % Angle of ray where each robot has final position
koopman_file = 'dubin_learned_koopman.mat';         % File containing learned Koopman model
r_margin = 0.14;                                    % Minimum distance between robot center points                
alpha = 1;                                          % CBF strengthening term
obs = [0;0];                                        % Center of obstacle
r_circ = 0.2;                                       % Radius of obstacle (- r_margin)

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
options = qpOASES_options('printLevel',0);              % Solver options for supervisory controller
affine_dynamics = @(x) dubin(x);                        % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                        % State is defined as x = [X,Y,v,theta], u = [a,r]
dt = r.time_step;
u_lim = [-0.1 0.1; - 2*pi/3, 2*pi/3];
barrier_func_collision = @(x_1, x_2) collision_avoidance_vec(x_1,x_2,r_margin);                   
supervisory_controller = @(x,u0,agent_ind) koopman_qp_cbf_multiagent_vec(x, u0, agent_ind, N_max, affine_dynamics, barrier_func_collision, alpha, N, func_dict, K_pows, C, options, u_lim);

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
    if exist('sc1','var')
        delete(sc1);
    end
    sc1 = scatter(final(1,:),final(2,:),40,1:N,'filled');
    
    while(~init_checker(x, final))
        
        x = r.get_poses();
        if exist('sc2','var')
            delete(sc2);
        end
        sc2 = scatter(x(1,:),x(2,:),40,1:N);
        
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
        
        for i = 1:N
            u_barrier = supervisory_controller(x_mod,u_0(:,i),i);
            vd_est = v_est(i) + u_barrier(1)*dt;
        
            dxu(:,i) = [vd_est; u_barrier(2)];
        end
        
        r.set_velocities(1:N, dxu);
        r.step();   
        v_prev = dxu(1,:);
        x_prev = x;
        
    end    
end

% We can call this function to debug our experiment!  Fix all the errors
% before submitting to maximize the chance that your experiment runs
% successfully.
r.debug();