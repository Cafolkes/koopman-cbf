% Initializing the agents to random positions with barrier certificates. 
% This script shows how to initialize robots to a particular point.
% Sean Wilson
% 07/2019
close all; clc; clear all;

n = 3;
N = 1;

% Generate intial and final positions:
radius = 0.8;
start_angles = linspace(0,2*pi-(2*pi/N),N);
end_angles = wrapTo2Pi(start_angles + pi);

initial_positions = zeros(n,N);
initial_positions(1,:) = radius*cos(start_angles);
initial_positions(2,:) = radius*sin(start_angles);
initial_positions(3,:) = start_angles-pi;

final_positions = zeros(n,N);
final_positions(1,:) = radius*cos(end_angles);
final_positions(2,:) = radius*sin(end_angles);

% Construct Koopman CBF supervisory controller:
load('dubin_learned_koopman.mat')
options = qpOASES_options('printLevel',0);          % Solver options for supervisory controller
affine_dynamics = @(x) dubin(x);                    % System dynamics, returns [f,g] with x_dot = f(x) + g(x)u
                                                    % State is defined as x = [X,Y,v,theta], u = [a,r]
r_margin = 0.05;                                    % Minimum distance between robot center points                
barrier_func_collision = @(x_1, x_2) collision_avoidance(x_1,x_2,r_margin);                   
obs = [0;0];
r_circ = 0.1;
alpha = 1;
barrier_func = @(x) round_obs(x,obs,r_circ);             % Barrier function
                                                    
%supervisory_controller = @(x,u0) koopman_qp_cbf_multiagent(x, u0, N_max, affine_dynamics, barrier_func_collision, N, func_dict, K_pows, C, options);
supervisory_controller = @(x,u0) koopman_qp_cbf_static(x, u0, N_max, affine_dynamics, barrier_func, alpha, func_dict, K_pows, C, options);

%initial_positions = generate_initial_conditions(N, 'Spacing', 0.5);
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', initial_positions);

% Create a barrier certificate so that the robots don't collide
%si_barrier_certificate = create_si_barrier_certificate();
si_to_uni_dynamics = create_si_to_uni_dynamics();

% We'll make the rotation error huge so that the initialization checker
% doesn't care about it
args = {'PositionError', 0.02, 'RotationError', 50};
init_checker = create_is_initialized(args{:});
controller = create_si_position_controller();

% Get initial location data for while loop condition.
x=r.get_poses();
r.step();

n_passes = 2;
v_prev = zeros(1,N);
x_prev = x;
dt = r.time_step;
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
        dxi = controller(x(1:2, :), final(1:2, :)); % Legacy controller (greedy position residual)
        dxu = si_to_uni_dynamics(dxi, x);
        
        v_x_est = (x(1,:)-x_prev(1,:))/dt;
        v_y_est = (x(2,:)-x_prev(2,:))/dt;
        v_est = sqrt(v_x_est.^2+v_y_est.^2);
        x_mod = [x(1,:); x(2,:); v_est; x(3,:)];
        
        ad_est = (dxu(1,:) - v_est)/dt;
        u_0 = [ad_est; dxu(2,:)];
        
        %for i = 1:N
            % TODO: See how controller was solved in previous paper..
            %u_barrier = supervisory_controller(x_mod, u_0(:,i));
        %    u_barrier = supervisory_controller(x_mod(:,i),u_0(:,i));
        %    vd_est = v_est(i) + u_barrier(1)*dt;
        %    dxi(:,i) = [vd_est; u_barrier(2)];
        %end
        
        u_barrier = supervisory_controller(x_mod,u_0);
        vd_est = v_est + u_barrier(1)*dt;
        
        if dxu ~= [vd_est; u_barrier(2)]
            fprintf('CBF active! diff =%f\n', norm(dxu-[vd_est; u_barrier(2)]))
        end
        dxu = [vd_est; u_barrier(2)];
   
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

