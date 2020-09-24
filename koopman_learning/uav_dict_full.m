
x = sym('x',[17,1],'real');
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
x17 = x(17);

p = [x1;x2;x3];
q = [x4;x5;x6;x7];
v = [x8;x9;x10];
w = [x11;x12;x13];
Omega = [x14;x15;x16;x17];

D = [1;p;q;v;w;Omega];


% Polynomials up to 2nd degree of quaternion and angular velocity terms:
qw = [q;w];
for i=1:length(qw)
    for j=i:length(qw)
        D = [D; qw(i)*qw(j)];
    end
end

%Add Omega squared times body z-direction
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

%Add Omega squared terms and omega squared times body in z-direction
q_sq = q(1)^2+q(2)^2+q(3)^2+q(4)^2;
for i = 1:length(Omega)
    D = [D; Omega(i)^2];
    D = [D; zBodyInWorld/q_sq*Omega(i)^2];
end

% Polynomials of quaternion states divided by ||q||^2
%q_sq = q(1)^2+q(2)^2+q(3)^2+q(4)^2;
%for i=1:length(q)
%    for j=i:length(q)
%        D = [D; q(i)*q(j)/q_sq];
%    end
%end

% Polynomials up to 3rd degree of quaternion and angular velocity states:
% qw = [q;w];
% for i=1:length(qw)
%     for j=i:length(qw)
%         for k=j:length(qw)
%             D = [D; qw(i)*qw(j)*qw(k)];
%         end
%     end
% end

D = [D;
    Omega(1)*w(3)
    atan2(2*q(1)*q(2) + 2*q(3)*q(4), - 2*q(2)^2 - 2*q(3)^2 + 1);
    atan2(2*q(1)*q(2) + 2*q(3)*q(4), - 2*q(2)^2 - 2*q(3)^2 + 1);
    (v(2)*(2*q(3)^2 + 2*q(4)^2 - 1))/((2*q(3)^2 + 2*q(4)^2 - 1)^2 + (2*q(1)*q(4) + 2*q(2)*q(3))^2)^(1/2);
    (v(1)*(2*q(1)*q(4) + 2*q(2)*q(3)))/((2*q(3)^2 + 2*q(4)^2 - 1)^2 + (2*q(1)*q(4) + 2*q(2)*q(3))^2)^(1/2);
    asin(2*q(1)*q(3) - 2*q(2)*q(4));
    (v(1)*(2*q(3)^2 + 2*q(4)^2 - 1))/((2*q(3)^2 + 2*q(4)^2 - 1)^2 + (2*q(1)*q(4) + 2*q(2)*q(3))^2)^(1/2);
    (v(2)*(2*q(1)*q(4) + 2*q(2)*q(3)))/((2*q(3)^2 + 2*q(4)^2 - 1)^2 + (2*q(1)*q(4) + 2*q(2)*q(3))^2)^(1/2);
    1/(2*q(2)^2 + 2*q(3)^2 - 1)];


J = jacobian(D,x);

F = [D,J];
n_lift = length(D);
C=[zeros(17,1) eye(17) zeros(17,n_lift-1-17)];

matlabFunction(D,J,C,'file','koopman_learning/uav_D_full.m');