function u=MIQP_MPC_v3(x0,des,N)
global Ts vm am rm
if nargin<3
N=20;
end
if size(des,1)==1
    des=des';
end
psi=x0(4);
vx=x0(3);
desired_psi=atan_2pi(des(1)-x0(1),des(2)-x0(2));
if abs(desired_psi-psi)>pi
    if psi>0
        desired_psi=desired_psi+2*pi;
    else
        desired_psi=desired_psi-2*pi;
    end
end
des(4)=desired_psi-psi;
Ac=[0 0 cos(psi) -vx*sin(psi);0 0 sin(psi) vx*cos(psi);0 0 0 0;0 0 0 0];
Bc=[0 0;0 0;1 0;0 1];
[Ad,Bd]=c2d(Ac,Bc,Ts);
x0(4)=0;
n=size(Ad,1);
m=size(Bd,2);
[Lx,Lu]=mpc_mat(Ad,Bd,N);
t=Ts*(1:N)';
alpha=0.1;
Q_final=diag([1 1 0 0.1]);
Q_process=diag([1 1 0 1]);
Q_process=kron(diag(exp(alpha*(-N:-1))),Q_process);
R=0.01*diag(repmat([1 1],1,N));

A_v=[0 0 1 0];
b_v=vm;
A_mpc_v=[];
b_mpc_v=[];
for i=1:N
    A_mpc_v=[A_mpc_v;A_v*Lu((i-1)*n+1:i*n,:)];
    b_mpc_v=[b_mpc_v;b_v-A_v*Lx((i-1)*n+1:i*n,:)*x0];

end
A_u=[];
b_u=[];
for i=1:N
    A_u=blkdiag(A_u,[eye(2);-eye(2)]);
    b_u=[b_u;[am;rm;am;rm]];
end
Q_des_final=Lu((N-1)*n+1:N*n,:)'*Q_final*Lu((N-1)*n+1:N*n,:);
Q_des_process=Lu'*Q_process*Lu;

f_des_final=(Lx((N-1)*n+1:N*n,:)*x0-des)'*Q_final*Lu((N-1)*n+1:N*n,:);
f_des_process=(x0'*Lx'-repmat(des,N,1)')*Q_process*Lu ;

Q_des=Q_des_final+Q_des_process*1/N;
f_des=f_des_final+f_des_process*1/N;

tol=1e-10;
Aineq=[A_mpc_v;A_u];
bineq=[b_mpc_v;b_u];
% options = optimoptions('quadprog', 'Algorithm','interior-point-convex','Display','off');
% [u,~,flag]=quadprog(Q_des+R,f_des,Aineq,bineq);
options = qpOASES_options('printLevel',0);
%     options.printLevel = PL_LOW;
[sol,~,exitflag] =qpOASES( Q_des+R,f_des',Aineq,[],[],[],bineq,options);
u=sol(1:2);


    