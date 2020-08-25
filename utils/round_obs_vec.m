function b = round_obs_vec(x,center,r)
    b = diag((x(:,1:2)-center)*(x(:,1:2)-center)')-r^2;
end