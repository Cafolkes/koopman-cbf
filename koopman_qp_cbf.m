function u = koopman_qp_cbf(x, u0, N, system_dynamics, func_dict, K_pows, C, options)
    global Ts alpha obs r

    [d,J] = func_dict(x);
    xx = zeros(N,4);
    QQ = zeros(N*4,4);
    tt = Ts*(1:N);
    for j=1:N
        xx(j,:)=(C*K_pows{j}*d)';
        QQ(4*(j-1)+1:4*j,:)=C*K_pows{j}*J;
    end
    Aineq = [];
    bineq = [];
    [f,g] = system_dynamics(x);
    for j = 1:length(tt)
        b = round_obs(xx(j,:)',obs,r);
        if b<1
            h = 1e-4;
            db = zeros(4,1);
            db(1) = (round_obs(xx(j,:)'+[h;0;0;0],obs,r)-b)/h;
            db(2) = (round_obs(xx(j,:)'+[0;h;0;0],obs,r)-b)/h;
            db(3) = (round_obs(xx(j,:)'+[0;0;h;0],obs,r)-b)/h;
            db(4) = (round_obs(xx(j,:)'+[0;0;0;h],obs,r)-b)/h;
            Aineq = [Aineq;-db'*QQ(4*(j-1)+1:4*j,:)*g];
            bineq = [bineq;alpha*b+db'*QQ(4*(j-1)+1:4*j,:)*f];
        end
    end
    if isempty(Aineq)
        u = u0;
    else
        [u,~,~] =qpOASES(eye(2),-u0,Aineq,[],[],[],bineq,options);
    end
end


function [b,db] = round_obs(x,center,r)
    b = (x(1:2)-center)'*(x(1:2)-center)-r^2;
    db = [2*(x(1:2)-center);zeros(2,1)];
end