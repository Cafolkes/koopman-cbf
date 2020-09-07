function u = koopman_qp_cbf_multiagent_obstacle_vec_mod(x, u0, agent_ind, N, system_dynamics, backup_dynamics, barrier_func_obstacle, barrier_func_collision, alpha, n_agents, func_dict, CK_pows, options, u_lim,n,m)
 assert(size(x,2)==n_agents,'Number of agents misspesified');
    
    xx = zeros(n_agents,N,n);
    h = 1e-4;
    
    for i = 1 : n_agents
        [d,J] = func_dict(x(:,i));
        for k=1:N
            xx(i,k,:)=(CK_pows{k}*d)';
            QQ{i,k} = CK_pows{k}*J;
        end
    end
    
    Aineq = [];
    bineq = [];
    
    f_cl = backup_dynamics(x(:,agent_ind));
    [f,g] = system_dynamics(x(:,agent_ind));
    x_i = reshape(xx(agent_ind,:,:),N,n);
    
    % Static obstacle avoidance:
    b_obs = barrier_func_obstacle(x_i);

    inds = find(b_obs < 1);
    x_i_red = x_i(inds,:);
    qq = blkdiag(QQ{agent_ind,inds});
    b_obs_red = b_obs(inds);
    N_red = size(inds,1);    

    if N_red > 0
        db_obs = zeros(N_red,n);
        for k = 1 : n
            x_pert = zeros(1,n);
            x_pert(k) = h;
            db_obs(:,k) = (barrier_func_obstacle(x_i_red+x_pert)-b_obs_red)/h;
        end

        db_obs_cell = mat2cell(db_obs,ones(1,N_red));
        db_obs_blkdiag = blkdiag(db_obs_cell{:});
        f_ext = repmat((f-f_cl),N_red,1);
        g_ext = repmat(g,N_red,1);
        Aineq = [Aineq;-db_obs_blkdiag*qq*g_ext -ones(N_red,1)];
        bineq = [bineq;alpha*b_obs_red+db_obs_blkdiag*qq*f_ext];
    end
    
    % Agent-agent collision avoidance:
    for j = 1 : n_agents
        if agent_ind == j
            continue
        end

        x_j = reshape(xx(j,:,:),N,n);
        b_coll = barrier_func_collision(x_i, x_j);

        inds = find(b_coll < 1);
        x_i_red = x_i(inds,:);
        x_j_red = x_j(inds,:);
        qq = blkdiag(QQ{agent_ind,inds});
        b_coll_red = b_coll(inds);
        N_red = size(inds,1);

        if N_red > 0
            db_coll = zeros(N_red,n);
            for k = 1 : n
                x_pert = zeros(1,n);
                x_pert(k) = h;
                db_coll(:,k) = (barrier_func_collision(x_i_red+x_pert,x_j_red)-b_coll_red)/h;
            end
            
            db_coll_cell = mat2cell(db_coll,ones(1,N_red));
            db_coll_blkdiag = blkdiag(db_coll_cell{:});
            f_ext = repmat((f-f_cl),N_red,1);
            g_ext = repmat(g,N_red,1);
            Aineq = [Aineq;-db_coll_blkdiag*qq*g_ext -ones(N_red,1)];
            bineq = [bineq;alpha*b_coll_red+db_coll_blkdiag*qq*f_ext];
        end
    end
    
    nonzero_inds = find(all(Aineq(:,1:m),2));
    Aineq = Aineq(nonzero_inds,:);
    bineq = bineq(nonzero_inds);
    
    if isempty(Aineq)
        u = u0;
    else
        H = diag([ones(1,m) 0]);
        [res,~,~] = quadprog(H,[-u0;1e5],Aineq,bineq,[],[],[u_lim(:,1);0],[u_lim(:,2);inf],[u0;0],options);
        %[u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
        u = res(1:m);
        %if res(m+1) > 1e-2
            %disp(strcat('infeasible, agent',num2str(agent_ind)))
        %end
    end
end