function b = collision_avoidance_vec(x_1,x_2,r_margin)
    b = diag((x_1(:,1:2)-x_2(:,1:2))*(x_1(:,1:2)-x_2(:,1:2))')-r_margin^2;
end