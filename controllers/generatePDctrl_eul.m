%% Generate PD-controller
clc;
generatePDcontroller()

function generatePDcontroller()
% define states
x = sym('x',[3,1],'real');
v = sym('v',[3,1],'real');
eul = sym('eul',[3,1],'real');
w = sym('w',[3,1],'real');
Omega = sym('Omega',[4,1],'real');
X = [x;eul;v;w;Omega];
Kpxy = sym('Kpxy',[1,1],'real');
Kpz = sym('Kpz',[1,1],'real');
KpVxy = sym('KpVxy',[1,1],'real');
KpVz = sym('KpVz',[1,1],'real');
KpAtt = sym('KpAtt',[1,1],'real');
KdAtt = sym('KdAtt',[1,1],'real');
KpOmegaz = sym('KpOmegaz',[1,1],'real');
hoverT = sym('hoverT',[1,1],'real');
M = [Kpxy;Kpz;KpVxy;KpVz;KpAtt;KdAtt;KpOmegaz;hoverT];

x_d = sym('x_d',[3,1],'real');
v_d = sym('v_d',[3,1],'real');
eul_d = sym('q_d',[3,1],'real');
w_d = sym('w_d',[3,1],'real');
X_d = [x_d;eul_d;v_d;w_d];

rotm = [cos(eul(3))*cos(eul(2)), cos(eul(3))*sin(eul(1))*sin(eul(2)) - cos(eul(1))*sin(eul(3)), sin(eul(1))*sin(eul(3)) + cos(eul(1))*cos(eul(3))*sin(eul(2));
   cos(eul(2))*sin(eul(3)), cos(eul(1))*cos(eul(3)) + sin(eul(1))*sin(eul(3))*sin(eul(2)), cos(eul(1))*sin(eul(3))*sin(eul(2)) - cos(eul(3))*sin(eul(1));
    -sin(eul(2)),          cos(eul(2))*sin(eul(1)),                              cos(eul(1))*cos(eul(2))];

zBodyInWorld = rotm*[0;0;1];

% velocity
vWorldNoYaw = [cos(-eul(3)) -sin(-eul(3)) 0; sin(-eul(3)) cos(-eul(3)) 0; 0 0 1]*v;

xError = Kpxy*(x_d(1) - x(1));
yError = Kpxy*(x_d(2) - x(2));
zError = x_d(3) - x(3);

vxError = KpVxy*(v_d(1)-vWorldNoYaw(1));
vyError = KpVxy*(v_d(2)-vWorldNoYaw(2));
vzError = v_d(3)-vWorldNoYaw(3);

rollError = eul_d(1)-vyError-yError-eul(1);
pitchError = eul_d(2)+vxError+xError-eul(2);
rollErrorDot = w_d(1)-w(1);
pitchErrorDot = w_d(2)-w(2);
yawErrorDot = w_d(3)-w(3);

uz = Kpz*zError+KpVz*vzError+hoverT/zBodyInWorld(3);
uroll = KpAtt*rollError+KdAtt*rollErrorDot;
upitch = KpAtt*pitchError+KdAtt*pitchErrorDot;
uyaw = KpOmegaz*yawErrorDot;

U(1) = uz-uroll-upitch-uyaw;
U(2) = uz-uroll+upitch+uyaw;
U(3) = uz+uroll+upitch-uyaw;
U(4) = uz+uroll-upitch+uyaw;
U = transpose(U);
matlabFunction(U,'Vars',{X,X_d,M},'File','pdU_eul');
end