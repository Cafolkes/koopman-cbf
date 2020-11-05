function plot_training_fit(X, K_pows, C, func_dict, error_bound)
    global Ts x_bdry
    
    %Calculate fit:
    for i = 1 : length(X)
        x0 = X{i}(1,:);
        [z0,~] = func_dict(x0);
        x_hat = [x0'];
        for j = 1 : size(X{i},1)-1
            x_hat = [x_hat C*K_pows{j}*z0];
        end
        X_hat{i} = x_hat';
    end
        
    fig = figure(2);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    
    subplot(1,2,1)
    ind_x = 3;
    ind_y = 4;
    hold on
    for i = 1 : length(X_hat)
        X{i}(:,4) = wrapTo2Pi(X{i}(:,4));
        scatter(X{i}(:,ind_x),X{i}(:,ind_y))
    end
    xlabel('Velocity ($x_3$)');
    ylabel('Angle ($x_4$)');
    xlim([0 x_bdry(ind_x,2)]);
    ylim([x_bdry(ind_y,1) x_bdry(ind_y,2)]);
    title('Training data (states)')
    
    
    err_bnd = error_bound(X{1}(1,:)');
    bound = [0];
    for k = 1 : length(X{1})-1
        bound = [bound err_bnd{k}];
    end
    
    subplot(1,2,2)    
    hold on
    tt = 0 : Ts : Ts*(length(X{i})-1);
    plot(tt, bound, '--r', 'lineWidth',2)
    for i = 1 : length(X_hat)
        tt = 0 : Ts : Ts*(size(X{i},1)-1);
        diff = X{i}-X_hat{i};
        diff(:,4) = angdiff(X{i}(:,4),X_hat{i}(:,4));
        plot(tt, vecnorm(diff,2,2))
    end
    xlabel('Time (sec)');
    ylabel('Prediction error $||x-\hat{x}||$');
    title('Training error')
    legend('Error bound')
    
    saveas(fig,'figures/training_fit.png')
    
end