function [u, int_time] = qp_cbf_obs(x, u0, N, system_dynamics, backup_dynamics, barrier_func_obs, alpha, sensitivity_dynamics, options,u_lim,n,m)
    global Ts
    t0 = posixtime(datetime('now'));
    w0 = [x; reshape(eye(n),n*n,1)];
    if N>= 1
        [~, w] = ode45(sensitivity_dynamics, [0:Ts:Ts*N], w0);
    else
        w = w0;
    end
    xx = w(1:N,1:n);
    QQ = zeros(N*n,n);
    for i = 1 : N
        QQ((i-1)*n+1:i*n,:) = reshape(w(i,n+1:end),n,n);
    end
    int_time = posixtime(datetime('now')) - t0; 
    
    Aineq = [];
    bineq = [];
    
    f_cl = backup_dynamics(0, x);
    [f,g] = system_dynamics(x);
    
    Aineq = zeros(N,m);
    bineq = zeros(N,1);
    for k = 1:N
        x_1 = reshape(xx(k,:),n,1);

        b = barrier_func_obs(x_1);
        qq = QQ(n*(k-1)+1:n*k,:);

        h = 1e-4;
        db = zeros(n,1);
        for l = 1 : n
            x_pert = zeros(n,1);
            x_pert(l) = h;
            db(l) = (barrier_func_obs(x_1+x_pert)-b)/h;
        end
        Aineq(k,:) = -db'*qq*g;
        bineq(k) = alpha*b+db'*qq*(f-f_cl);
    end
    
    if isempty(Aineq)
        u = u0;
    else
        nonzero_inds = find(all(Aineq(:,1:m),2));
        Aineq = Aineq(nonzero_inds,:);
        bineq = bineq(nonzero_inds);
         
        Aineq = [Aineq -ones(size(Aineq,1),1)];
        H = diag([ones(1,m) 0]);
        [res,~,~] = quadprog(H,[-u0;1e6],Aineq,bineq,[],[],[u_lim(:,1);0],[u_lim(:,2);inf],[u0;0],options);
        
        u = res(1:m);
    end
end