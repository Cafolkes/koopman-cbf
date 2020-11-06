function J_cl = get_jacobian_cl(dynamics_cl, n, m)
    x = sym('x',[n,1],'real');
    f_cl = dynamics_cl(x);
    J_sym = jacobian(f_cl, x);
    J_cl = matlabFunction(J_sym, 'Vars', {x});
end