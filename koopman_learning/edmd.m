function [K, C] = edmd(X, func_dict)
    [Z, Z_p] = lift_data(X, func_dict);
    n_lift = size(Z,1);

    cvx_begin
    cvx_solver gurobi
    variables K(n_lift,n_lift)
    minimize(norm(Z_p - K*Z, Inf))
    cvx_end

    [~,~,C] = func_dict(X{1}(1,:));
end

function [Z, Z_p] = lift_data(X, func_dict)
Z = [];
Z_p = [];
for i = 1:size(X,2)
    for j = 2:size(X{i},1)
        [z,~] = func_dict(X{i}(j-1,:));
        [z_p,~] = func_dict(X{i}(j,:));
        Z = [Z z];
        Z_p = [Z_p z_p];
    end
end
end