function err_bnd = koopman_error_bound(x,X,L_f, L_phi,e_max,tt,K_pows,C,func_dict, bound_type)
    %Bound type: 0 -> local error, 1 -> global error, old, 2 -> global
    %error, modified
    global Ts;
    
    x(4) = wrapTo2Pi(x(4)); % Shift all angles to be between [0,2pi]
    x_hat = find_closest_x(x,X);
    x_hat(4) = wrapTo2Pi(x_hat(4)); % Shift all angles to be between [0,2pi]
    p=2;
    d = func_dict(x);
    d_hat = func_dict(x_hat);
    
    CA_norm = 0;
    diff_x = x-x_hat;
    diff_x(4) = angdiff(x(4),x_hat(4));
    
    err_bnd = {};
    for k = 1:length(tt)-1
        CA_norm = CA_norm + norm(C*K_pows{1},p)^(k); %Note: shifted because error for time=0 not calculated
        switch bound_type
            case 0
                err_bnd{k} = norm(e_max,p)*CA_norm;
            case 1
                err_bnd{k} = norm(C*K_pows{k}*(d-d_hat),p) + norm(e_max,p)*CA_norm + norm(diff_x,p);
            case 2
                A_hat = K_pows{1}-eye(size(K_pows{1},1));
                eps_max = L_f*Ts*norm(diff_x,p) + e_max + norm(C*A_hat*(d-d_hat),p);
                err_bnd{k} = norm(eps_max,p)*CA_norm;
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