function err_bnd = koopman_error_bound_mu(mu,L,e_max,tt,K_pows,C,func_dict, bound_type)
    
    global Ts;
    
    p=2;    
    CA_norm = 0;
    
    err_bnd = {};
    for k = 1:length(tt)-1
        CA_norm = CA_norm + norm(C*K_pows{k},p); %Note: shifted because error for time=0 not calculated
        switch bound_type
            case 0
                err_bnd{k} = norm(e_max,p)*CA_norm;
            case 1
                err_bnd{k} = norm(C*K_pows{k},p)*L*mu + norm(e_max,p)*CA_norm + mu;
        end
                
    end 
end

function x_hat = find_closest_x(x,X)
    n_traj = length(X);
    dist = [];
    min_res = 1e6;
    traj = 1;
    samp = 1;
    for i = 1 : n_traj
        X_mod = X{i};   
        X_mod(:,4) = wrapTo2Pi(X_mod(:,4));
        %dist = [dist; vecnorm(X_mod-x',2,2)];
        dist = vecnorm(X_mod-x',2,2);
        [min_tmp, argmin_tmp] = min(dist);
        if min_tmp < min_res
            traj = i;
            samp = argmin_tmp;
            min_res = min_tmp;
        end
    end
    x_hat = X{traj}(samp,:)';
end