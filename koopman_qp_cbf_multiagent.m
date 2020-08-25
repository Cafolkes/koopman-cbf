function u = koopman_qp_cbf_multiagent(x, u0, N, system_dynamics, barrier_func_collision, alpha, n_agents, func_dict, K_pows, C, options)
    assert(size(x,2)==n_agents,'Number of agents misspesified');
    n=4;
    
    xx = zeros(n_agents,N,n);
    h = 1e-4;
    
    for i = 1 : n_agents
        [d,J] = func_dict(x(:,i));
        for k=1:N
            xx(i,k,:)=(C*K_pows{k}*d)';
            %QQ{i}(n*(k-1)+1:n*k,:)=C*K_pows{k}*J;
            QQ{i,k} = C*K_pows{k}*J;
        end
        %QQ_blkdiag{i} = blkdiag(QQ{i,:});
    end
    %db_blkdiag_mult = 
    
    Aineq = [];
    bineq = [];
    
    for i = 1 : n_agents
        [f,g] = system_dynamics(x(:,i));
        %qq = QQ_blkdiag{i};
        x_1 = reshape(xx(i,:,:),N,n);
        for j = 1 : n_agents
            if i == j
                continue
            end
            
            x_2 = reshape(xx(j,:,:),N,n);
            
            %for k = 1:N
                %x_1 = reshape(xx(i,k,:),n,1);
                %x_2 = reshape(xx(j,k,:),n,1);
                
            b = barrier_func_collision(x_1, x_2);
            
            inds = find(b < 1);
            x_1_red = x_1(inds,:);
            x_2_red = x_2(inds,:);
            qq = blkdiag(QQ{i,inds});
            b_red = b(inds);
            N_red = size(inds,1);
            
            if N_red > 0
                db = zeros(N_red,n);
                db(:,1) = (barrier_func_collision(x_1_red+[h 0 0 0],x_2_red)-b_red)/h;
                db(:,2) = (barrier_func_collision(x_1_red+[0 h 0 0],x_2_red)-b_red)/h;
                db(:,3) = (barrier_func_collision(x_1_red+[0 0 h 0],x_2_red)-b_red)/h;
                db(:,4) = (barrier_func_collision(x_1_red+[0 0 0 h],x_2_red)-b_red)/h;

                db_ext = repmat(db,1,N_red);
                f_ext = repmat(f,N_red,1);
                g_ext = repmat(g,N_red,1);
                Aineq = [Aineq;-db_ext*qq*g_ext];
                bineq = [bineq;alpha*b_red+db_ext*qq*f_ext];
            end
             %   end
            %end
        end
    end
    
    if isempty(Aineq)
        u = u0;
    else
        [u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
    end
end