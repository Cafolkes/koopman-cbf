function [Z, Z_p] = lift_data(X, func_dict, center_data, first_obs_one)
Z = [];
Z_p = [];
n_lift = length(func_dict(zeros(size(X{1},1))));
for i = 1:size(X,2)
    x = X{i};
    z_end = func_dict(x(end,:));
    z_end(5) = z_end(5) + 1e-1*randn();
    for j = 2:size(X{i},1)
        z = func_dict(x(j-1,:));
        z_p = func_dict(x(j,:));
        
        if center_data == true
            %z_end(5) = 0;
            if first_obs_one ==true
                z_end(1) = 0;
            end
            z = z - z_end;
            z_p = z_p - z_end;
        end
        
        Z = [Z z];
        Z_p = [Z_p z_p];
    end
end
end