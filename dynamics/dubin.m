function [f,g] = dubin(x)

f = [x(3)*cos(x(4));x(3)*sin(x(4));0;0];
g = [zeros(2,2);eye(2)];

end