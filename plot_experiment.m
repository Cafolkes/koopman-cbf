function plot_experiment(x_rec, u_rec, u0_rec, func_dict, K_pows, C)

    global Ts am T_exp obs r
    
    for i=1:round(T_exp/Ts)
        figure(3)
        clf
        hold on
        draw_circle(obs(1),obs(2),r);
        x = x_rec(i,:);
        plot(x(1),x(2),'rO');
        quiver(x(1),x(2),0.2*x(3)*cos(x(4)),0.2*x(3)*sin(x(4)));
        [d,~] = func_dict(x);
        N = ceil(x(3)/am/Ts);
        xx = zeros(N,4);
        for j=1:N
            xx(j,:)=(C*K_pows{j}*d)';
        end
        for j = 1:N
            plot(xx(j,1),xx(j,2),'b*');
        end
        axis([-2,15,-5,5])
        axis equal
        drawnow
    end
end