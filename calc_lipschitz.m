function L = calc_lipschitz(n,m,affine_dynamics, con1)
    global x_bdry
    x = sym('x',[n,1]);
    [f,g] = affine_dynamics(x);
    dynamics_cl = f+g*con1(x);
    J_cl = jacobian(dynamics_cl,x);
    
    L_row = zeros(n,1);
    lb = x_bdry(:,1);
    ub = x_bdry(:,2);
    
    for i = 1:n
        J_fun = matlabFunction(J_cl(i,:),'Vars',x);
        switch n
            case 2
                J_fun = @(x) J_fun(x(1),x(2));
            case 3
                J_fun = @(x) J_fun(x(1),x(2),x(3));
            case 4
                J_fun = @(x) J_fun(x(1),x(2),x(3),x(4));
            case 5
                J_fun = @(x) J_fun(x(1),x(2),x(3),x(4),x(5));
            case 6
                J_fun = @(x) J_fun(x(1),x(2),x(3),x(4),x(6));
            otherwise
                disp('Warning, higher state dimension than 6 not implemented')
                return
            
        obj = @(x) -norm(J_fun(x),2)
        
        x_tmp = fmincon(obj,rand(4,1),[],[],[],[],lb,ub);
        
        L_row(i) = -obj(x_tmp);
    end
    L = max(L_row);
end