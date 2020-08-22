function plot_training_fit(X, K_pows, C, func_dict, error_bound)
    global Ts
    
    %Calculate fit:
    for i = 1 : length(X)
        x0 = X{i}(1,:);
        [z0,~] = func_dict(x0);
        z = [x0'];
        for j = 1 : length(X{i})-1
            z = [z C*K_pows{j}*z0];
        end
        X_hat{i} = z;
    end
        
    fig = figure(1);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    
    subplot(1,2,1)
    ind_x = 3;
    ind_y = 4;
    hold on
    for i = 1 : length(X_hat)
        scatter(X{i}(ind_x,:),X{i}(ind_y,:))
    end
    xlabel('Velocity ($x_3$)');
    ylabel('Angle ($x_4$)');
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
        tt = 0 : Ts : Ts*(length(X{i})-1);
        plot(tt, vecnorm(X{i}-X_hat{i}',2,2))
    end
    xlabel('Time (sec)');
    ylabel('Prediction error $||x-\hat{x}||$');
    title('Training error')
    legend('Error bound')
    
    saveas(fig,'figures/training_fit.png')
    
    
    
end