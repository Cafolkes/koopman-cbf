function u = koopman_qp_cbf_multiagent(x, u0, N, system_dynamics, barrier_func_collision, alpha, n_agents, func_dict, K_pows, C, options)
    assert(size(x,2)==n_agents,'Number of agents misspesified');
    n=4;
    
    xx = zeros(n_agents,N,n);
    
    for i = 1 : n_agents
        [d,J] = func_dict(x(:,i));
        for k=1:N
            xx(i,k,:)=(C*K_pows{k}*d)';
            QQ{i}(n*(k-1)+1:n*k,:)=C*K_pows{k}*J;
        end
    end
    
    Aineq = [];
    bineq = [];
    
    for i = 1 : n_agents
        [f,g] = system_dynamics(x(:,i));
        for j = 1 : n_agents
            if i == j
                continue
            end
            
            for k = 1:N
                x_1 = reshape(xx(i,k,:),n,1);
                x_2 = reshape(xx(j,k,:),n,1);
                
                b = barrier_func_collision(x_1, x_2);
                if b<1
                    qq = QQ{i}(n*(k-1)+1:n*k,:);
                    
                    h = 1e-4;
                    db = zeros(n,1);
                    db(1) = (barrier_func_collision(x_1+[h;0;0;0],x_2)-b)/h;
                    db(2) = (barrier_func_collision(x_1+[0;h;0;0],x_2)-b)/h;
                    db(3) = (barrier_func_collision(x_1+[0;0;h;0],x_2)-b)/h;
                    db(4) = (barrier_func_collision(x_1+[0;0;0;h],x_2)-b)/h;
                    Aineq = [Aineq;-db'*qq*g];
                    bineq = [bineq;alpha*b+db'*qq*f];
                end
            end
        end
    end
    
    if isempty(Aineq)
        u = u0;
    else
        [u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
    end
end