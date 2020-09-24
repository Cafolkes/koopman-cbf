function b = collision_avoidance_3d_vec(x_1,x_2,r_margin)
    %b = diag((x_1(:,1:3)-x_2(:,1:3))*(x_1(:,1:3)-x_2(:,1:3))')-r_margin^2;
    b = vecnorm(x_1(:,1:3)-x_2(:,1:3),2,2).^2-r_margin^2;
end