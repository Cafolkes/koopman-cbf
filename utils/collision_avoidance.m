function [b,db] = collision_avoidance(x_1,x_2,r_margin)
    b = (x_1(1:2)-x_2(1:2))'*(x_1(1:2)-x_2(1:2))-r_margin^2;
    db = [2*(x_1(1:2)-x_2(1:2));zeros(2,1);-2*(x_1(1:2)-x_2(1:2));zeros(2,1);];
end