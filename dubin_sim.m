function xdot = dubin_sim(x,u)
    global vm
    [f,g] = dubin(x);
    xdot = f+g*u;
    if x(3)<-vm
        xdot(3)=max(u(1),-x(3)-vm);
    elseif x(3)>vm
        xdot(3)=min(u(1),vm-x(3));
    end
end