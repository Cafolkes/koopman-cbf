x = sym('x',[4,1]);
X = x(1);
Y = x(2);
v = x(3);
psi = x(4);
D = [1;X;Y;v;psi;v^2;cos(psi);sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
     v^3*cos(psi);v^3*sin(psi)];
 
J = jacobian(D,x);

F = [D,J];
C=[zeros(4,1) eye(4) zeros(4,9)];
matlabFunction(D,J,C,'file','dubin_D.m');
