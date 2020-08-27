function [D,J,C] = dubin_D(x1,x2,x3,x4)
%DUBIN_D
%    [D,J,C] = DUBIN_D(X1,X2,X3,X4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    27-Aug-2020 17:54:15

t2 = cos(x4);
t3 = sin(x4);
t4 = x3.^2;
t5 = t3.*x3;
t6 = t2.*x3;
t7 = t3.*t4;
t8 = t2.*t4;
t9 = t3.*t4.*x3;
t10 = t2.*t4.*x3;
D = [1.0;x1;x2;x3;x4;t4;t2;t3;t6;t5;t8;t7;t10;t9];
if nargout > 1
    J = reshape([0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,x3.*2.0,0.0,0.0,t2,t3,t2.*x3.*2.0,t3.*x3.*2.0,t2.*t4.*3.0,t3.*t4.*3.0,0.0,0.0,0.0,0.0,1.0,0.0,-t3,t2,-t5,t6,-t7,t8,-t9,t10],[14,4]);
end
if nargout > 2
    C = reshape([0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,14]);
end
