function e_max = calc_max_residual(X, func_dict, K, C)
    e_vec = [];
    for i = 1:length(X)
        x = X{i};
        x_pred=[];
        for j = 1:size(x,1)-1
            z = func_dict(x(j,:));
            x_pred = [x_pred C*K*z];
        end
        res = x(2:end,:) - x_pred';
        e_vec = [e_vec vecnorm(res,2,2)];
    end
    e_max = max(max(e_vec));
end
    