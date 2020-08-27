function xdot = sim_cl_dynamics(x, u, affine_dynamics)
    [f,g] = affine_dynamics(x);
    xdot = f + g*u;
end