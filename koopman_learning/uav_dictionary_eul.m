
x = sym('x',[16,1],'real');
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);
x8 = x(8);
x9 = x(9);
x10 = x(10);
x11 = x(11);
x12 = x(12);
x13 = x(13);
x14 = x(14);
x15 = x(15);
x16 = x(16);

p = [x1;x2;x3];
eul = [x4;x5;x6];
v = [x7;x8;x9];
w = [x10;x11;x12];
Omega = [x13;x14;x15;x16];

D = [1;p;eul;v;w;Omega];

% Ang vel terms:
D = [D;
    w(3)*cos(eul(1))*tan(eul(2));
    w(2)*sin(eul(1))*tan(eul(2));
    w(2)*cos(eul(1));
    w(3)*sin(eul(1));
    (w(3)*cos(eul(1)))/cos(eul(2));
    (w(2)*sin(eul(1)))/cos(eul(2))];

% Linear acceleration terms:
rotm = [cos(eul(3))*cos(eul(2)), cos(eul(3))*sin(eul(1))*sin(eul(2)) - cos(eul(1))*sin(eul(3)), sin(eul(1))*sin(eul(3)) + cos(eul(1))*cos(eul(3))*sin(eul(2));
   cos(eul(2))*sin(eul(3)), cos(eul(1))*cos(eul(3)) + sin(eul(1))*sin(eul(3))*sin(eul(2)), cos(eul(1))*sin(eul(3))*sin(eul(2)) - cos(eul(3))*sin(eul(1));
    -sin(eul(2)),          cos(eul(2))*sin(eul(1)),                              cos(eul(1))*cos(eul(2))];
D = [D;
    rotm*[0;0;Omega(1)^2];
    rotm*[0;0;Omega(2)^2];
    rotm*[0;0;Omega(3)^2];
    rotm*[0;0;Omega(4)^2]];
    
% Angular acceleration terms:
D = [D;
    Omega(1)^2;
    Omega(2)^2;
    Omega(3)^2;
    Omega(4)^2;
    w(2)*w(3);
    w(1)*w(3);
    w(1)*w(2)];

D = [D;
    v(1)*cos(eul(3));
    v(2)*sin(eul(3));
    v(2)*cos(eul(3));
    v(1)*sin(eul(3));
    1/(cos(eul(1))*cos(eul(2)))];
    
J = jacobian(D,x);

F = [D,J];
n_lift = length(D);
C=[zeros(16,1) eye(16) zeros(16,n_lift-1-16)];

matlabFunction(D,J,C,'file','koopman_learning/uav_D_eul.m');