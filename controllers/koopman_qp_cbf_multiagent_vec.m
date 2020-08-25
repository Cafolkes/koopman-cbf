function u = koopman_qp_cbf_multiagent_vec(x, u0, agent_ind, N, system_dynamics, barrier_func_collision, alpha, n_agents, func_dict, K_pows, C, options, u_lim)
    assert(size(x,2)==n_agents,'Number of agents misspesified');
    n=4;
    
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

        inds = find(b < 0.5);
        x_1_red = x_i(inds,:);
        x_2_red = x_j(inds,:);
        qq = blkdiag(QQ{agent_ind,inds});
        b_red = b(inds);
        N_red = size(inds,1);

        if N_red > 0
            db = zeros(N_red,n);
            db(:,1) = (barrier_func_collision(x_1_red+[h 0 0 0],x_2_red)-b_red)/h;
            db(:,2) = (barrier_func_collision(x_1_red+[0 h 0 0],x_2_red)-b_red)/h;
            db(:,3) = (barrier_func_collision(x_1_red+[0 0 h 0],x_2_red)-b_red)/h;
            db(:,4) = (barrier_func_collision(x_1_red+[0 0 0 h],x_2_red)-b_red)/h;

            db_cell = mat2cell(db,ones(1,N_red));
            db_blkdiag = blkdiag(db_cell{:});
            f_ext = repmat(f,N_red,1);
            g_ext = repmat(g,N_red,1);
            Aineq = [Aineq;-db_blkdiag*qq*g_ext];
            bineq = [bineq;alpha*b_red+db_blkdiag*qq*f_ext];
        end
    end
    
    if nargin > 12
       Aineq = [Aineq;-eye(2);eye(2)]; % [other constraints; lower lim, upper lim]
       bineq = [bineq;-u_lim(:,1);u_lim(:,2)];
    end
    
    if isempty(Aineq)
        u = u0;
    else
        [u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
    end
end