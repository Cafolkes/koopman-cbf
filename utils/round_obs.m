function [b,db] = round_obs(x,center,r)
    b = (x(1:2)-center)'*(x(1:2)-center)-r^2;
    db = [2*(x(1:2)-center);zeros(2,1)];
end