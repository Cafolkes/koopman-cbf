function dubin_dictionary(add_bias)
    x = sym('x',[4,1]);
    X = x(1);
    Y = x(2);
    v = x(3);
    psi = x(4);
    %
    %D = [X;Y;v;psi;v^2;cos(psi)-1;sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
    %     v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
    %D = [1;X;Y;v;psi;v^2;v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
    %     v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
    if add_bias == true
        %D = [1;X;Y;v;psi;v^2;cos(psi);sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
        %     v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
        D = [1; X;Y;v;psi;v^2;v^3;v^4;v^5;cos(psi);sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
             v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
        C = [zeros(4,1) eye(4) zeros(4,16)];
        %C = [zeros(4,1) eye(4) zeros(4,13)];
        fname = 'koopman_learning/dubin_D.m';
    else
        D = [X;Y;v;psi;v^2;cos(psi);sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
             v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
        D = [X;Y;v;psi;v^2;v^3;v^4;v^5;cos(psi);sin(psi);v*cos(psi);v*sin(psi);v^2*cos(psi);v^2*sin(psi);...
             v^3*cos(psi);v^3*sin(psi); v^4*cos(psi);v^4*sin(psi); v^5*cos(psi);v^5*sin(psi)];
        C = [eye(4) zeros(4,16)];
        %C = [eye(4) zeros(4,13)];
        fname = 'koopman_learning/dubin_D_nb.m';
    end
    J = jacobian(D,x);

    %F = [D,J];
    %C=[zeros(4,1) eye(4) zeros(4,9)];
    %C=[zeros(4,1) eye(4) zeros(4,13)];
    
    %C=[eye(4) zeros(4,16)];
    %C=[zeros(4,1) eye(4) zeros(4,11)];
    matlabFunction(D,J,C,'file',fname);
end
