
x = sym('x',[13,1],'real');
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

p = [x1;x2;x3];
q = [x4;x5;x6;x7];
v = [x8;x9;x10];
w = [x11;x12;x13];

D = [1;p;q;v;w];


% Polynomials up to 2nd degree excluding position terms (dynamics are
% position invariant):
qvw = [q;v;w];
for i=1:length(qvw)
    for j=i:length(qvw)
        D = [D; qvw(i)*qvw(j)];
    end
end

% Polynomials of quaternion states divided by ||q||^2
q_sq = q(1)^2+q(2)^2+q(3)^2+q(4)^2;
for i=1:length(q)
    for j=i:length(q)
        D = [D; q(i)*q(j)/q_sq];
    end
end

% Polynomials up to 3rd degree of quaternion and angular velocity states:
qw = [q;w];
for i=1:length(qw)
    for j=i:length(qw)
        for k=j:length(qw)
            D = [D; qw(i)*qw(j)*qw(k)];
        end
    end
end


J = jacobian(D,x);

F = [D,J];
n_lift = length(D);
C=[zeros(13,1) eye(13) zeros(13,n_lift-1-13)];

matlabFunction(D,J,C,'file','koopman_learning/uav_D.m');