%% Compute gradients of dynamics
addpath('dynamics')
generateLinearizedDynamics;
%%
function generateLinearizedDynamics()
clear all;
con = quad1_constants();

x = sym('x',[3,1],'real');
v = sym('v',[3,1],'real');
q = sym('q',[4,1],'real');
w = sym('w',[3,1],'real');
Omega = sym('Omega',[4,1],'real');

X = [x;q;v;w;Omega];
KpVxy = sym('KpVxy',[1,1],'real');
KpVz = sym('KpVz',[1,1],'real');
KpAtt = sym('KpAtt',[1,1],'real');
KdAtt = sym('KdAtt',[1,1],'real');
KpOmegaz = sym('KpOmegaz',[1,1],'real');
hoverT = sym('hoverT',[1,1],'real');
M = [KpVxy;KpVz;KpAtt;KdAtt;KpOmegaz;hoverT];
V=backupU(X,M);
d = 0;

[dr, ddr, dxi_q, dw, dOmega] = full_eul(x,v,q,w,Omega,V*14.8,d,con.m,con.g,con.J_bod,con.D,con.J_rot, con.J_prop, con.K_v, con.R, con.k_f, con.k_t,true);
Xdot = [dr;dxi_q;ddr;dw;dOmega];

A_lin = jacobian(Xdot,X);
matlabFunction(A_lin,'Vars',{X,M},'File','matAlin');

Omega = zeros(4,1);
[dr, ddr, dxi_q, dw, ~] = full_eul(x,v,q,w,Omega,V*14.8,d,con.m,con.g,con.J_bod,con.D,con.J_rot, con.J_prop, con.K_v, con.R, con.k_f, con.k_t,true);
Xdot_red = [dr;dxi_q;ddr;dw];
X_red = [x;q;v;w];

A_lin_red = jacobian(Xdot_red,X_red);
matlabFunction(A_lin_red,'Vars',{X_red,M},'File','matAred');
end