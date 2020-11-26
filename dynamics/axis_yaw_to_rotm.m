%% axis_yaw_to_rotm: convert b3c and yaw to rotation matrix
function [R, dR, d2R] = axis_yaw_to_rotm(b3c, yaw, db3c, dyaw, d2b3c, d2yaw)

	if nargin < 6 
		d2yaw = 0;
		d2b3c = [0;0;0];
	end

	if nargin < 4
		dyaw = 0;
		db3c = [0;0;0];
	end

	Rz = [cos(yaw) -sin(yaw) 0; 
		  sin(yaw)  cos(yaw) 0;
		  0         0        1];

	Rz1 = [-sin(yaw) -cos(yaw) 0; 
 		   cos(yaw)  -sin(yaw) 0;
 		   0         0        0];

    Rz2 = [-cos(yaw) sin(yaw) 0;
            -sin(yaw) -cos(yaw) 0;
            0         0         0];

    dRz = Rz1 * dyaw;
    d2Rz = Rz2 * dyaw^2 + Rz1 * d2yaw;

	x = Rz'*b3c;
	dx = dRz'*b3c + Rz'*db3c;
	d2x = d2Rz'*b3c + 2*dRz'*db3c + Rz'*d2b3c;

	c = sqrt(x(1)^2 + x(3)^2);
	dc = (x(1) * dx(1) + x(3) * dx(3)) / c;
	d2c = (x(1) * d2x(1) + dx(1)^2 + x(3) * d2x(3) + dx(3)^2)/c ...
	      - (x(1) * dx(1) + x(3) * dx(3))^2/c^3;

	if x(3) < 0
		% make sure cos(theta) >= 0  (pitch angle goes -pi/2 to pi/2)
		c = -c;
		dc = -dc;
		d2c = -d2c;
	end

	ci = 1/c;
	dci = -(1/c^2)*dc;
	d2ci = (2/c^3)*dc^2 - (1/c^2)*d2c;

	x1x2 = x(1)*x(2);
	dx1x2 = dx(1)*x(2) + x(1)*dx(2);
	d2x1x2 = d2x(1)*x(2) + 2*dx(1)*dx(2) + x(1)*d2x(2);

	x2x3 = x(2)*x(3);
	dx2x3 = dx(2)*x(3) + x(2)*dx(3);
	d2x2x3 = d2x(2)*x(3) + 2*dx(2)*dx(3) + x(2)*d2x(3);

	RyRx = [x(3)*ci  -x1x2*ci  x(1);
			0         c        x(2);
			-x(1)*ci -x2x3*ci  x(3)];

	dRyRx = zeros(3,3);
	dRyRx(1,1) = dx(3)*ci + x(3)*dci;
	dRyRx(1,2) = -dx1x2*ci - x1x2*dci;
	dRyRx(1,3) = dx(1);
	dRyRx(2,1) = 0;
	dRyRx(2,2) = dc;
	dRyRx(2,3) = dx(2);
	dRyRx(3,1) = -dx(1)*ci - x(1)*dci;
	dRyRx(3,2) = -dx2x3*ci - x2x3*dci;
	dRyRx(3,3) = dx(3);

	d2RyRx = zeros(3,3);
	d2RyRx(1,1) = d2x(3)*ci + 2*dx(3)*dci + x(3)*d2ci;
	d2RyRx(1,2) = -d2x1x2*ci - 2*dx1x2*dci - x1x2*d2ci;
	d2RyRx(1,3) = d2x(1);
	d2RyRx(2,1) = 0;
	d2RyRx(2,2) = d2c;
	d2RyRx(2,3) = d2x(2);
	d2RyRx(3,1) = -d2x(1)*ci - 2*dx(1)*dci - x(1)*d2ci;
	d2RyRx(3,2) = -d2x2x3*ci - 2*dx2x3*dci - x2x3*d2ci;
	d2RyRx(3,3) = d2x(3);

	% final rotation matrix
	R = Rz * RyRx;
	dR = dRz * RyRx + Rz * dRyRx;
	d2R = d2Rz * RyRx + 2 * dRz * dRyRx + Rz * d2RyRx;