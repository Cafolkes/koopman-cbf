function [D,J,C] = dubin_D(x1,x2,x3,x4)
%DUBIN_D
%    [D,J,C] = DUBIN_D(X1,X2,X3,X4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Nov-2020 22:13:30

t2 = cos(x4);
t3 = sin(x4);
t4 = x3.^2;
t5 = x3.^3;
t7 = x3.^5;
t6 = t4.^2;
t8 = t2.*x3;
t9 = t3.*x3;
t10 = t2.*t4;
t11 = t2.*t5;
t13 = t2.*t7;
t14 = t3.*t4;
t15 = t3.*t5;
t17 = t3.*t7;
t12 = t2.*t6;
t16 = t3.*t6;
D = [1.0;x1;x2;x3;x4;t4;t5;t6;t7;t2;t3;t8;t9;t10;t14;t11;t15;t12;t16;t13;t17];
if nargout > 1
    J = reshape([0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,x3.*2.0,t4.*3.0,t5.*4.0,t6.*5.0,0.0,0.0,t2,t3,t8.*2.0,t9.*2.0,t10.*3.0,t14.*3.0,t11.*4.0,t15.*4.0,t12.*5.0,t16.*5.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,-t3,t2,-t9,t8,-t14,t10,-t15,t11,-t16,t12,-t17,t13],[21,4]);
end
if nargout > 2
    C = reshape([0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,21]);
end
