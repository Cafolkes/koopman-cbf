function u = koopman_qp_cbf_multiagent_vec(x, u0, agent_ind, N, system_dynamics, barrier_func_collision, alpha, n_agents, func_dict, K_pows, C, options, u_lim, n, m)
    assert(size(x,2)==n_agents,'Number of agents misspesified');
    xx = zeros(n_agents,N,n);
    h = 1e-4;
    
    for i = 1 : n_agents
        [d,J] = func_dict(x(:,i));
        for k=1:N
            xx(i,k,:)=(C*K_pows{k}*d)';
            QQ{i,k} = C*K_pows{k}*J;
        end
    end
    
    Aineq = [];
    bineq = [];
    
    [f,g] = system_dynamics(x(:,agent_ind));
    x_i = reshape(xx(agent_ind,:,:),N,n);
    for j = 1 : n_agents
        if agent_ind == j
            continue
        end

        x_j = reshape(xx(j,:,:),N,n);
        b = barrier_func_collision(x_i, x_j);

        inds = find(b < 1);
        x_1_red = x_i(inds,:);
        x_2_red = x_j(inds,:);
        qq = blkdiag(QQ{agent_ind,inds});
        b_red = b(inds);
        N_red = size(inds,1);

        if N_red > 0
            db = zeros(N_red,n);
            for k = 1 : n
                x_pert = zeros(1,n);
                x_pert(k) = h;
                db(:,k) = (barrier_func_collision(x_1_red+x_pert,x_2_red)-b_red)/h;
            end
            db_cell = mat2cell(db,ones(1,N_red));
            db_blkdiag = blkdiag(db_cell{:});
            f_ext = repmat(f,N_red,1);
            g_ext = repmat(g,N_red,1);
            Aineq = [Aineq;-db_blkdiag*qq*g_ext -ones(N_red,1)];
            %Aineq = [Aineq;-db_blkdiag*qq*g_ext];
            bineq = [bineq;alpha*b_red+db_blkdiag*qq*f_ext];
        end
    end
    
    if isempty(Aineq)
        u = u0;
    else
        H = diag([1 1 1 1 0]);
        [res,~,~] = quadprog(H,[-u0;1e4],Aineq,bineq,[],[],[u_lim(:,1);0],[u_lim(:,2);inf],[u0;0],options);
        %[u,~,flag] = qpOASES(eye(m),-u0,Aineq,u_lim(:,1),u_lim(:,2),[],bineq,options);
        %if flag ~= 0
        %    disp('Infeasible')
        %end
        u = res(1:4);
    end
end