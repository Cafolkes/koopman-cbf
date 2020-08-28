function [Z, Z_p] = lift_data(X, func_dict, nom_matrix)
Z = [];
Z_p = [];
for i = 1:size(X,2)
    for j = 2:size(X{i},1)
        z = func_dict(X{i}(j-1,1:13));
        z_p = func_dict(X{i}(j,1:13));
        
        if nargin > 2 %Subtract nominal model from z_p
            z_p = z_p - nom_matrix*z;
            %z_p = z_p - nom_matrix*z;
        end
        
        Z = [Z z];
        Z_p = [Z_p z_p];
    end
end
end