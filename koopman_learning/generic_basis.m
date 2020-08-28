function generic_basis(n,polyorder,useone)

    x = sym('x',[n,1],'real');

    D = [];

    % Constant term
    if useone == true
        D = [D; 1];
    end

    % Poly order 1
    if polyorder >= 1
        D = [D; x];
    end

    % Poly order 2
    if polyorder>=2
        for i=1:n
            for j=i:n
                D = [D; x(i)*x(j)];
            end
        end
    end

    % Poly order 3
    if polyorder>=3
        for i=1:n
            for j=i:n
                for k=j:n
                    D = [D; x(i)*x(j)*x(k)];
                end
            end
        end
    end

    if(polyorder>=4)
        disp('Warning: polyorder higher than 3 not implemented');
    end

    J = jacobian(D,x);

    F = [D,J];
    n_lift = length(D);
    if useone == true
        C=[zeros(n,1) eye(n) zeros(n,n_lift-1-n)];
    else
        C = [eye(n) zeros(n,n_lift-n)];
    end

    matlabFunction(D,J,C,'file',['koopman_learning/poly_D_' num2str(n) '_' num2str(polyorder) '.m']);
end