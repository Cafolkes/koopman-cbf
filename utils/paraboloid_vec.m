function h = paraboloid_vec(x,scale_fac,offset)
    h = x(:,3)-offset-scale_fac*(x(:,1).^2+x(:,2).^2);
end
