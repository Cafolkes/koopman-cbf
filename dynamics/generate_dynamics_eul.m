%% Compute gradients of dynamics
addpath('dynamics')
generateMatlabFunctions()

%% Compute backup controller and gradients
clc;
generateBackup()

%%
function generateMatlabFunctions()
clear all;
con = quad1_constants();

x = sym('x',[3,1],'real');
v = sym('v',[3,1],'real');
eul = sym('eul',[3,1],'real');
w = sym('w',[3,1],'real');
Omega = sym('Omega',[4,1],'real');
V=sym('V',[4,1],'real');
d = 0;
[dr, ddr, deul, dw, dOmega] = full_eul(x,v,eul,w,Omega,V*14.8,d,con.m,con.g,con.J_bod,con.D,con.J_rot, con.J_prop, con.K_v, con.R, con.k_f, con.k_t,false);
X = [x;eul;v;w;Omega];
Xdot = [dr;deul;ddr;dw;dOmega];
f = Xdot - jacobian(Xdot,V)*V;
g = jacobian(Xdot,V);
matlabFunction(f,'Vars',{X},'File','matf_eul');
matlabFunction(g,'Vars',{X},'File','matg_eul');
DFcl = jacobian(Xdot,X);
DF = jacobian(f,X);
DG = sym('DG', [length(X),length(V),length(X)],'real');

for i = 1:length(V)
    DG(:,i,:) = jacobian(g(:,i),X);
end
matlabFunction(DFcl,'Vars',{X,V},'File','matDFcl_eul');
matlabFunction(DF,'Vars',{X},'File','matDF_eul');
matlabFunction(DG,'Vars',{X},'File','matDG_eul');
end

function generateBackup()
% define states
eul_sequence = "XYZ";
x = sym('x',[3,1],'real');
v = sym('v',[3,1],'real');
eul = sym('eul',[3,1],'real');
w = sym('w',[3,1],'real');
Omega = sym('Omega',[4,1],'real');
X = [x;eul;v;w;Omega];
KpVxy = sym('KpVxy',[1,1],'real');
KpVz = sym('KpVz',[1,1],'real');
KpAtt = sym('KpAtt',[1,1],'real');
KdAtt = sym('KdAtt',[1,1],'real');
KpOmegaz = sym('KpOmegaz',[1,1],'real');
hoverT = sym('hoverT',[1,1],'real');
M = [KpVxy;KpVz;KpAtt;KdAtt;KpOmegaz;hoverT];

rotm = [cos(eul(3))*cos(eul(2)), cos(eul(3))*sin(eul(1))*sin(eul(2)) - cos(eul(1))*sin(eul(3)), sin(eul(1))*sin(eul(3)) + cos(eul(1))*cos(eul(3))*sin(eul(2));
   cos(eul(2))*sin(eul(3)), cos(eul(1))*cos(eul(3)) + sin(eul(1))*sin(eul(3))*sin(eul(2)), cos(eul(1))*sin(eul(3))*sin(eul(2)) - cos(eul(3))*sin(eul(1));
    -sin(eul(2)),          cos(eul(2))*sin(eul(1)),                              cos(eul(1))*cos(eul(2))];

zBodyInWorld = rotm*[0;0;1];

% velocity
vWorldNoYaw = [cos(-eul(3)) -sin(-eul(3)) 0; sin(-eul(3)) cos(-eul(3)) 0; 0 0 1]*v;

vDes = [0;0;0;0];
vxError = KpVxy*(vDes(1)-vWorldNoYaw(1));
vyError = KpVxy*(vDes(2)-vWorldNoYaw(2));
vzError = vDes(3)-vWorldNoYaw(3);

rollError = -vyError-eul(1);
pitchError = vxError-eul(2);
rollErrorDot = -w(1);
pitchErrorDot = -w(2);
yawErrorDot = -w(3);


uz = KpVz*vzError+hoverT/zBodyInWorld(3);
uroll = KpAtt*rollError+KdAtt*rollErrorDot;
upitch = KpAtt*pitchError+KdAtt*pitchErrorDot;
uyaw = KpOmegaz*yawErrorDot;

U(1) = uz-uroll-upitch-uyaw;
U(2) = uz-uroll+upitch+uyaw;
U(3) = uz+uroll+upitch-uyaw;
U(4) = uz+uroll-upitch+uyaw;
U = transpose(U);
DU = jacobian(U,X);
matlabFunction(U,'Vars',{X,M},'File','backupU_eul');
matlabFunction(DU,'Vars',{X,M},'File','backupDU_eul');
end