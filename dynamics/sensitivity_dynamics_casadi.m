function wdot = sensitivity_dynamics_casadi(w, J_cl, backup_dynamics, n)
% w = [x; q], q = flatten(Q)

W = reshape(w(n+1:end),n,n);

wdot = [backup_dynamics;
    reshape(J_cl*W, n*n,1)];
    