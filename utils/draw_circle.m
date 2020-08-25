function [] =draw_circle(X,Y,R,color,linewidth)
if nargin<4
    color='r';
    
end
if nargin<5
    linewidth=1;
end
theta=(1:360)*pi/180;
for i=1:size(theta,2)
    x(i)=X+R*cos(theta(i));
    y(i)=Y+R*sin(theta(i));
end
plot(x,y,color,'linewidth',linewidth);
fill(x,y,color);