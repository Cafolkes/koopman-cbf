function [f,g] = UAVDynamics_eul(X)
    f = matf_eul(X);
    g = matg_eul(X);
end