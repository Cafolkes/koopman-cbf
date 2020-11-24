function [K, obj_vals, lambda_tuned] = edmd(Z, Z_p, solver, first_obs_one, lambda_1, lambda_2, normalize, tune_fit, n_folds)
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
                        K(i,2:end) = k(:,inf.IndexMinMSE);
                        obj_vals(i) = inf.MSE(inf.IndexMinMSE);
                        lambda_tuned(i) = inf.LambdaMinMSE;
                        K(i,1) = inf.Intercept(inf.IndexMinMSE);
                    else
                        [K(i,2:end), inf] = lasso(x, y,'Lambda',lambda_1(i)); 
                        obj_vals(i) = inf.MSE;
                        K(i,1) = inf.Intercept; 
                    end
                else
                    x = Z';
                    col_normalizer = vecnorm(x);
                    x = x./col_normalizer;
                    y = Z_p(i,:);
                    if tune_fit == true
                        [k, inf] = lasso(x, y,'CV',n_folds,'Intercept',false,'Options',statset('UseParallel',true)); 
                        K(i,:) = k(:,inf.IndexMinMSE);
                        obj_vals(i) = inf.MSE(inf.IndexMinMSE);
                        lambda_tuned(i) = inf.LambdaMinMSE;
                        %K(i,1) = inf.Intercept(inf.IndexMinMSE);
                    else
                        %[K(i,:), inf] = lasso(x, y,'Lambda',1e-3,'Intercept',false);
                        [K(i,:), inf] = lasso(x, y,'Lambda',1e-3);
                        obj_vals(i) = inf.MSE;
                    end
                    K = K./col_normalizer;
                end
                fprintf('Learned %i out of %i predictors\n', i, n_pred);
            end
            
        elseif strcmp(solver,'gurobi')
            x = Z;
            y = Z_p;
            if normalize == true
                col_normalizer = std(x')';
                col_normalizer(find(col_normalizer<=1e-6)) = 1;
                x = x./col_normalizer;
            end
            
            cvx_begin
            cvx_solver gurobi
            variables K(n_pred,n_lift)
            obj = 0;
            for i = 1 : n_pred
                obj = obj + norm(K(i,:)*x-y(i,:),2) + lambda_1*norm(K(i,:),1) + lambda_2*norm(K(i,:),2);
            end
            minimize(obj)
            cvx_end
            
            if normalize == true
                K = K./col_normalizer';
            end
        end
    end
end