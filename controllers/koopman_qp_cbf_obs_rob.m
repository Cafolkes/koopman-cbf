function [u, int_time] = koopman_qp_cbf_obs_rob(x, u0, N, system_dynamics, backup_dynamics, barrier_func_obs, alpha, func_dict, CK_pows, options,u_lim,n,m, error_bound)
    t0 = posixtime(datetime('now'));
    [d,J] = func_dict(x);
    xx = CK_pows*d;
    QQ = CK_pows*J;
    int_time = posixtime(datetime('now')) - t0; 
    
    f_cl = backup_dynamics(x);
    [f,g] = system_dynamics(x);
    
    Aineq = zeros(2*N,m);
    bineq = zeros(2*N,1);
    err_bnd = cell2mat(error_bound(reshape(xx(1:n),n,1)));
    
    for k = 1:N
        x_1 = reshape(xx((k-1)*n+1:k*n),n,1);
        x_1 = midrad(x_1, err_bnd(k));

        b = barrier_func_obs(x_1);
        qq = QQ(n*(k-1)+1:n*k,:);

        h = intval(1e-4);
        db = intval(zeros(n,1));
        
        for l = 1 : n
            x_pert = intval(zeros(n,1));
            x_pert(l) = h;
            db(l) = (barrier_func_obs(x_1+x_pert)-b)/h;
        end
        
        Aineq(k,:) = -inf(db)'*qq*g;
        bineq(k) = alpha*inf(b)+inf(db)'*qq*(f-f_cl);
        Aineq(k+N,:) = -sup(db)'*qq*g;
        bineq(k+N) = alpha*sup(b)+sup(db)'*qq*(f-f_cl);
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