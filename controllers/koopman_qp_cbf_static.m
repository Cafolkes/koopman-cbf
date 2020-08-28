function u = koopman_qp_cbf_static(x, u0, N, system_dynamics, barrier_func, alpha, func_dict, K_pows, C, options, u_lim)
    [d,J] = func_dict(x);
    xx = zeros(N,4);
    QQ = zeros(N*4,4);
    
    for j=1:N
        xx(j,:)=(C*K_pows{j}*d)';
        QQ(4*(j-1)+1:4*j,:)=C*K_pows{j}*J;
    end
    Aineq = [];
    bineq = [];
    [f,g] = system_dynamics(x);
    for j = 1:N
        b = barrier_func(xx(j,:)');
        if b<1
            h = 1e-4;
            db = zeros(4,1);
            db(1) = (barrier_func(xx(j,:)'+[h;0;0;0])-b)/h;
            db(2) = (barrier_func(xx(j,:)'+[0;h;0;0])-b)/h;
            db(3) = (barrier_func(xx(j,:)'+[0;0;h;0])-b)/h;
            db(4) = (barrier_func(xx(j,:)'+[0;0;0;h])-b)/h;
            Aineq = [Aineq;-db'*QQ(4*(j-1)+1:4*j,:)*g -1];
            bineq = [bineq;alpha*b+db'*QQ(4*(j-1)+1:4*j,:)*f];
        end
    end
    
    if isempty(Aineq)
        u = u0;
    else
        H = diag([1 1 0]);
        [res,~,~] = quadprog(H,[-u0;1e3],Aineq,bineq,[],[],[u_lim(:,1);0],[u_lim(:,2);inf],[u0;0],options);
        %[u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
        %[u,~,~] =qpOASES(eye(2),-u0,Aineq(:,1:2),[],[],[],bineq);
        u = res(1:2);
    end
end