function x0 = generate_initial_state_uav(use_quat)
    global x_bdry;
    
    n_sample = 16; % Number of states to sample
    x_sample = x_bdry(:,1)+(x_bdry(:,2)-x_bdry(:,1)).*rand(n_sample,1);
    if use_quat == true
        x0 = [x_sample(1:3); eul2quat(x_sample(4:6)')'; x_sample(7:end)];
    else
        x0 = x_sample;
    end
end