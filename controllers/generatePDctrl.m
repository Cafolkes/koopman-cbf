%% Generate PD-controller
clc;
generatePDcontroller()

function generatePDcontroller()
% define states
x = sym('x',[3,1],'real');
v = sym('v',[3,1],'real');
q = sym('q',[4,1],'real');
w = sym('w',[3,1],'real');
Omega = sym('Omega',[4,1],'real');
X = [x;q;v;w;Omega];
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
q_d = sym('q_d',[4,1],'real');
w_d = sym('w_d',[3,1],'real');
X_d = [x_d;q_d;v_d;w_d];

rotm(1,1) = 1-2*q(3)*q(3)-2*q(4)*q(4);
rotm(1,2) = 2*q(2)*q(3)-2*q(4)*q(1);
rotm(1,3) = 2*q(2)*q(4)+2*q(3)*q(1);
rotm(2,1) = 2*q(2)*q(3)+2*q(4)*q(1);
rotm(2,2) = 1-2*q(2)*q(2)-2*q(4)*q(4);
rotm(2,3) = 2*q(3)*q(4)-2*q(2)*q(1);
rotm(3,1) = 2*q(2)*q(4)-2*q(3)*q(1);
rotm(3,2) = 2*q(3)*q(4)+2*q(2)*q(1);
rotm(3,3) = 1-2*q(2)*q(2)-2*q(3)*q(3);

zBodyInWorld = rotm*[0;0;1];

% velocity
eul_tmp = 2.0 * (q(3) * q(3));
eulx = atan2(2*q(1)*q(2)+2*q(3)*q(4),(1-2*(q(2)*q(2))-eul_tmp));
euly = asin(2*(q(1)*q(3)-q(4)*q(2)));
eulz = atan2(2.0 * q(1) * q(4) + 2.0 * q(2) * q(3), (1.0 - eul_tmp) - 2.0 * (q(4) * q(4)));
vWorldNoYaw = [cos(-eulz) -sin(-eulz) 0; sin(-eulz) cos(-eulz) 0; 0 0 1]*v;

eul_tmp_d = 2.0 * (q_d(3) * q_d(3));
eulx_d = atan2(2*q_d(1)*q_d(2)+2*q_d(3)*q_d(4),(1-2*(q_d(2)*q_d(2))-eul_tmp_d));
euly_d = asin(2*(q_d(1)*q_d(3)-q_d(4)*q_d(2)));
%eulz_d = atan2(2.0 * q_d(1) * q_d(4) + 2.0 * q_d(2) * q_d(3), (1.0 - eul_tmp_d) - 2.0 * (q_d(4) * q_d(4)));

xError = Kpxy*(x_d(1) - x(1));
yError = Kpxy*(x_d(2) - x(2));
zError = x_d(3) - x(3);

vxError = KpVxy*(v_d(1)-vWorldNoYaw(1));
vyError = KpVxy*(v_d(2)-vWorldNoYaw(2));
vzError = v_d(3)-vWorldNoYaw(3);

rollError = eulx_d-vyError-yError-eulx;
pitchError = euly_d+vxError+xError-euly;
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
matlabFunction(U,'Vars',{X,X_d,M},'File','pdU');
end