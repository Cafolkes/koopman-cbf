%% get derivatives of desired z axis
function [b3c, db3c, d2b3c] = b3c_derivatives(b3c_pren, db3c_pren, d2b3c_pren)

	% calculate derivatives of (1/||.||)
	N = 1/norm(b3c_pren);
	dN = -(db3c_pren'*b3c_pren)/norm(b3c_pren)^3;
	d2N = -(d2b3c_pren'*b3c_pren + db3c_pren'*db3c_pren)/norm(b3c_pren)^3 ...
	      + 3*(db3c_pren'*b3c_pren)^2/norm(b3c_pren)^5;

	b3c = N * b3c_pren;
	db3c = dN * b3c_pren + N * db3c_pren;
	d2b3c = d2N * b3c_pren + 2*dN*db3c_pren + N*d2b3c_pren;
end
