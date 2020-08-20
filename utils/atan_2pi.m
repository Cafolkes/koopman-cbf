function theta=atan_2pi(x,y)
    if x>0
        theta=atan(y/x);
    elseif x==0
        theta=sign(y)*pi/2;
    else
        if y>0
            theta=atan(y/x)+pi;
        else
            theta=atan(y/x)-pi;
        end
end