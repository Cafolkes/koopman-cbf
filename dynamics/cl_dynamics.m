function f_cl = cl_dynamics(x,affine_dynamics, backup_controller)
   [f,g] = affine_dynamics(x);
   f_cl = f + g*backup_controller(x);
end