function f_cl = cl_dynamics(x,affine_dynamics, backup_controller, backup_controller_process)
   [f,g] = affine_dynamics(x);
   u = backup_controller(x);
   u = backup_controller_process(u);
   f_cl = f + g*u;
end