% Initializing the agents to random positions with barrier certificates. 
% This script shows how to initialize robots to a particular point.
% Sean Wilson
% 07/2019
close all; clc; 

% Experiment parameters:
n = 3;                                              % Number of states (unicycle model)
N = 5;                                              % Number of robots
n_passes = 2;                                       % Number of times to travel between waypoints
radius_waypoints = 0.8;                             % Radius of circle where waypoints are placed
start_angles = linspace(0,2*pi-(2*pi/N),N);         % Angle of ray where each robot is placed initially
end_angles = wrapTo2Pi(start_angles + pi);          % Angle of ray where each robot has final position
koopman_file = 'data/dubin_learned_koopman.mat';         % File containing learned Koopman model
r_margin = 0.08;                                    % Minimum distance between robot center points                
alpha = 1.25;                                          % CBF strengthening term
obs = [-0.3; 0.1];                                   % Center of obstacle
r_obs = 0.2;                                        % Radius of obstacle (- r_margin)

% Data storage:
file_name = 'data/collision_obstacle_avoidance.mat';     % File to save data matrices
for i = 1 : N
    x_data{i} = [];                                 % Store state of each robot
    backup_data{i} = [];                            % Store difference between legacy and supervisory controller (norm(u0-u))
    x_init_data{i} = [];                            % Store initial point
    x_final_data{i} = [];                           % Store final point
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
args = {'PositionError', 0.15, 'RotationError', 50};    
init_checker = create_is_initialized(args{:});          % Checks if desired position is reached
args = {'VelocityMagnitudeLimit', 0.15};    
controller = create_si_position_controller(args{:});           % Legacy controller (greedy position controller)

% Construct Koopman CBF supervisory controller:
load(koopman_file)
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));
options = optimoptions('quadprog','Display','none');    % Solver options for supervisory controller
affine_dynamics = @(x) dubin(x);                        % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                        % State is defined as x = [X,Y,v,theta], u = [a,r]
dt = r.time_step;
u_lim = [-0.1 0.1; - 2*pi/3, 2*pi/3];
barrier_func_collision = @(x_1, x_2) collision_avoidance(x_1,x_2,2*r_margin);                   
barrier_func_obstacle = @(x) round_obs(x,obs,r_obs);             
draw_circle(obs(1),obs(2),r_obs-r_margin);
func_dict = @(x) dubin_D(x(1),x(2),x(3),x(4));
rm = 3*pi/2;                                            % Maximum yaw rate
am = 0.1;                                               % Maximum acceleration
backup_controller = @(x) [-am*sign(x(3)); abs(x(3))*rm/vm];       % Backup controller (max brake and max turn rate)
backup_controller_process = @(u) u;
backup_dynamics = @(x) cl_dynamics(x,affine_dynamics, backup_controller, backup_controller_process);
supervisory_controller = @(x,u0,agent_ind) koopman_qp_cbf_multi_obs_coll(x, u0, agent_ind, N_max, affine_dynamics, backup_dynamics, barrier_func_collision, barrier_func_obstacle, alpha, N, func_dict, cell2mat(CK_pows'), options, u_lim, 4,2);
                                                                        
% Get initial location data for while loop condition.
x=r.get_poses();
r.step();

v_prev = zeros(1,N);
x_prev = x;
first_it = true;
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
    sc1 = scatter(final(1,:),final(2,:),100,1:N,'filled');
    
    while(~init_checker(x, final))
        
        x = r.get_poses();
        for i = 1:N
            if first_it == false
                delete(traj(i));
            end
            traj(i) = plot([x(1,i) final(1,i)],[x(2,i) final(2,i)],':k');
        end
        first_it = false;
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
            %if abs(dxu(1,i)) > pi
            %    disp(dxu(2,i))
            %end
            u_barrier = supervisory_controller(x_mod,u_0(:,i),i);
            vd_est = v_est(i) + u_barrier(1)*dt;
            
            vd_est = max(vd_est,0); % Ensure vd is in [0. 0.15] (as a result of the simplified state estimation)
            vd_est = min(vd_est,0.15);
            dxu(:,i) = [vd_est; u_barrier(2)];
        end
        
        % Threshold values to avoid actuator errors:
        max_frac_lin = 3/4;
        max_frac_ang = 1/4;
        dxu(1,:) = min(dxu(1,:),ones(1,N)*r.max_linear_velocity*max_frac_lin);
        dxu(1,:) = max(dxu(1,:),zeros(1,N));
        dxu(2,:) = min(dxu(2,:),ones(1,N)*r.max_angular_velocity*max_frac_ang);
        dxu(2,:) = max(dxu(2,:),-ones(1,N)*r.max_angular_velocity*max_frac_ang);
        
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