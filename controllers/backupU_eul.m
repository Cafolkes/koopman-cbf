function U = backupU_eul(in1,in2)
%BACKUPU_EUL
%    U = BACKUPU_EUL(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Aug-2020 19:25:55

KdAtt1 = in2(4,:);
KpAtt1 = in2(3,:);
KpOmegaz1 = in2(5,:);
KpVz1 = in2(2,:);
KpVxy1 = in2(1,:);
eul1 = in1(4,:);
eul2 = in1(5,:);
eul3 = in1(6,:);
hoverT1 = in2(6,:);
v1 = in1(7,:);
v2 = in1(8,:);
v3 = in1(9,:);
w1 = in1(10,:);
w2 = in1(11,:);
w3 = in1(12,:);
t2 = cos(eul1);
t3 = cos(eul2);
t4 = cos(eul3);
t5 = sin(eul3);
t6 = KpVz1.*v3;
t7 = KdAtt1.*w1;
t8 = KdAtt1.*w2;
t9 = KpOmegaz1.*w3;
t10 = t4.*v1;
t11 = t4.*v2;
t12 = t5.*v1;
t13 = t5.*v2;
t14 = -t6;
t15 = -t7;
t16 = -t8;
t17 = -t9;
t18 = 1.0./t2;
t19 = 1.0./t3;
t20 = -t12;
t21 = t10+t13;
t24 = hoverT1.*t18.*t19;
t22 = KpVxy1.*t21;
t23 = t11+t20;
t25 = eul2+t22;
t26 = KpVxy1.*t23;
t27 = -t26;
t28 = KpAtt1.*t25;
t29 = eul1+t27;
t30 = -t28;
t31 = KpAtt1.*t29;
t32 = -t31;
U = [t7+t8+t9+t14+t24+t28+t31;t7+t14+t16+t17+t24+t30+t31;t9+t14+t15+t16+t24+t30+t32;t8+t14+t15+t17+t24+t28+t32];