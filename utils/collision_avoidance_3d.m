function b = collision_avoidance_3d(x_1,x_2,r_margin)
    b = norm(x_1(1:3)-x_2(1:3))^2-r_margin^2;
end