function [K, obj_vals, lambda_tuned] = edmd(Z, Z_p, solver, first_obs_one, lambda, tune_fit, n_folds)
    n_lift = size(Z,1);
    n_pred = size(Z_p,1);
    obj_vals = zeros(n_pred,1);
    lambda_tuned = zeros(n_pred,1);
    
    if nargin <= 2
        cvx_begin
        cvx_solver gurobi
        variables K(n_lift,n_lift)
        minimize(norm(Z_p - K*Z, Inf))
        cvx_end
    else
        if strcmp(solver,'lasso')
            K = zeros(n_pred,n_lift);
            for i = 1 : n_pred
                if first_obs_one
                    x = Z(2:end,:)';
                    y = Z_p(i,:);
                    if tune_fit == true
                        [k, inf] = lasso(x, y,'CV',n_folds,'Options',statset('UseParallel',true)); 
                        K(i,1) = inf.Intercept(inf.IndexMinMSE);
                        K(i,2:end) = k(:,inf.IndexMinMSE);
                        obj_vals(i) = inf.MSE(inf.IndexMinMSE);
                        lambda_tuned(i) = inf.LambdaMinMSE;
                    else
                        [K(i,2:end), inf] = lasso(x, y,'Lambda',lambda(i)); 
                        K(i,1) = inf.Intercept;
                        obj_vals(i) = inf.MSE;
                    end
                else
                    x = Z(1:end,:)';
                    y = Z_p(i,:);
                    [K(i,:), inf] = lasso(x, y,'Lambda',1e-3,'Intercept',false);  
                    obj_vals(i) = inf.MSE;
                end
                fprintf('Learned %i out of %i predictors\n', i, n_pred);
            end
            
        elseif strcmp(solver,'gurobi')
            cvx_begin
            cvx_solver gurobi
            variables K(n_pred,n_lift)
            %obj = norm(Z_p - K*Z, Inf);
            obj = 0;
            for i = 1 : n_pred
                obj = obj + norm(K(i,:)*Z-Z_p(i,:),Inf) + 1e-3*norm(K(i,:),1);
            end
            minimize(obj)
            cvx_end
        end
    end
end