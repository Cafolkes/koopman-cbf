function u = koopman_qp_cbf_multiagent(x, u0, agent_ind, N, system_dynamics, backup_dynamics, barrier_func_collision, alpha, n_agents, func_dict, CK_pows, options,u_lim,n,m)
    assert(size(x,2)==n_agents,'Number of agents misspesified');
    
    xx = zeros(n_agents,N,n);
    
    for i = 1 : n_agents
        [d,J] = func_dict(x(:,i));
        for k=1:N
            xx(i,k,:)=(CK_pows{k}*d)';
            QQ{i}(n*(k-1)+1:n*k,:)=CK_pows{k}*J;
        end
    end
    
    Aineq = [];
    bineq = [];
    
    f_cl = backup_dynamics(x(:,agent_ind));
    [f,g] = system_dynamics(x(:,agent_ind));
    for j = 1 : n_agents
        if agent_ind == j
            continue
        end

        for k = 1:N
            x_1 = reshape(xx(agent_ind,k,:),n,1);
            x_2 = reshape(xx(j,k,:),n,1);

            b = barrier_func_collision(x_1, x_2);
            if b<1
                qq = QQ{agent_ind}(n*(k-1)+1:n*k,:);

                h = 1e-4;
                db = zeros(n,1);
                for l = 1 : n
                    x_pert = zeros(n,1);
                    x_pert(l) = h;
                    db(l) = (barrier_func_collision(x_1+x_pert,x_2)-b)/h;
                end
                Aineq = [Aineq;-db'*qq*g];
                bineq = [bineq;alpha*b+db'*qq*(f-f_cl)];
            end
        end
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
        %[u,~,~] =qpOASES(eye(m),-u0,Aineq,u_lim(:,1),u_lim(:,2),[],bineq,options);
    end
end