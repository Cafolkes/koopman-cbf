function err_bnd = koopman_error_bound(x,X,L,e_max,tt,K_pows,C,func_dict)
    global Ts
    
    x_hat = find_closest_x(x,X);
    p=2;
    [d,~] = func_dict(x);
    [d_hat,~] = func_dict(x_hat);
    
    CA_norm = 0;
    
    err_bnd = {};
    for k = 1:length(tt)-1
        if k >= 2
            CA_norm = CA_norm + norm(C*K_pows{1},p)^(k-1);
        end
        err_bnd{k} = norm(C*K_pows{k}*(d-d_hat),p) + norm(e_max,p)*CA_norm + exp(L*k*Ts)*norm(x-x_hat,p);
    end 
end

function x_hat = find_closest_x(x,X)
    n_traj = length(X);
    n_samp = length(X{1});
    dist = zeros(n_traj*n_samp,1);
    for i = 1 : n_traj
        dist((i-1)*n_samp+1:i*n_samp) = vecnorm(X{i}-x',2,2);
    end
    [~, argmin] = min(dist);
    traj = floor((argmin-1)/n_samp)+1;
    samp = argmin-((traj-1)*n_samp);
    x_hat = X{traj}(samp,:)';
end