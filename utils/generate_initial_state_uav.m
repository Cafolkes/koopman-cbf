function x0 = generate_initial_state_uav()
    global x_bdry;
    
    n_sample = 16; % Number of states to sample
    x_sample = x_bdry(:,1)+(x_bdry(:,2)-x_bdry(:,1)).*rand(n_sample,1);
    x0 = [x_sample(1:3); eul2quat(x_sample(4:6)')'; x_sample(7:end)];
    %x0 = [0; 0; 1; eul2quat(zeros(1,3))'; 1; 0; 0; zeros(3,1); 450*ones(4,1)]; %[p quat v w Omega]
end