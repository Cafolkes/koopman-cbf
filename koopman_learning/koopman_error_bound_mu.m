function err_bnd = koopman_error_bound_mu(mu,L_f, L_phi,e_max,tt,K_pows,C,func_dict, bound_type)
    %Bound type: 0 -> local error, 1 -> global error, old, 2 -> global
    %error, modified
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
                %err_bnd{k} = norm(C*K_pows{k}*mu,p) + norm(e_max,p)*CA_norm + norm(mu,p);
                err_bnd{k} = norm(C*K_pows{k}*mu,p) + norm(e_max,p)*CA_norm;
            case 2
                A_hat = K_pows{1}-eye(size(K_pows{1},1));
                %eps_max = L_f*Ts*mu + e_max + L_phi*norm(C*A_hat,p)*mu;
                eps_max = e_max + L_phi*norm(C*A_hat,p)*mu;
                err_bnd{k} = eps_max*CA_norm;
            case 3
                %err_bnd{k} = norm(C*K_pows{k}*mu,p) + norm(e_max,p)*CA_norm + norm(mu,p);
                err_bnd{k} = norm(C*(K_pows{k}-eye(size(K_pows{1},1)))*mu, p) + norm(e_max,p)*CA_norm;
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