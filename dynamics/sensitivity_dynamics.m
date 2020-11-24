function wdot = sensitivity_dynamics(w, J_cl, backup_dynamics, n)
% w = [x; q], q = flatten(Q)

W = reshape(w(n+1:end),n,n);

wdot = [backup_dynamics(w(1:n));
    reshape(J_cl(w(1:n))*W, n*n,1)];
    