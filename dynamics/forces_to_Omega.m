%% functionname: function description
function [Omega] = forces_to_Omega(fz, tau, D, k_f, k_t, Roc)

  if nargin == 6
    tau = tau + cross(Roc, [0;0;fz]);  % compensate for center of gravity
  end

  l = D/(2*sqrt(2));
  A = [  k_f    k_f    k_f    k_f;
      -l*k_f -l*k_f  l*k_f  l*k_f; 
      -l*k_f  l*k_f  l*k_f -l*k_f;
        -k_t    k_t   -k_t    k_t];
  
  Omega2 = A\[fz; tau];
  Omega = sign(Omega2) .* sqrt(abs(Omega2));
end