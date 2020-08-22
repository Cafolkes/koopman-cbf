function plot_test_fit(X_train,X_test, K_pows, C, func_dict, error_bound)
    global Ts
    
    %Calculate fit:
    for i = 1 : length(X_test)
        x0 = X_test{i}(1,:);
        [z0,~] = func_dict(x0);
        z = [x0'];
        for j = 1 : length(X_test{i})-1
            z = [z C*K_pows{j}*z0];
        end
        X_hat{i} = z;
    end
        
    fig = figure(2);
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    
    subplot(1,2,1)
    ind_x = 3;
    ind_y = 4;
    hold on
    for i = 1 : length(X_train)
        scatter(X_train{i}(ind_x,:),X_train{i}(ind_y,:),'r')
    end
    for i = 1 : length(X_test)
        scatter(X_test{i}(ind_x,:),X_test{i}(ind_y,:),'b')
    end
    legend('Training data', 'Test data');
    xlabel('Velocity ($x_3$)');
    ylabel('Angle ($x_4$)');
    title('Test data (states)');
    
    traj_bounds = [];
    for i = 1 : length(X_test)
        err_bnd = error_bound(X_test{i}(1,:)');
        bound = [0];
        for k = 1 : length(X_test{i})-1
            bound = [bound err_bnd{k}];
        end
        traj_bounds = [traj_bounds; bound];
    end
    [~,max_err_ind] = max(traj_bounds(:,2));
    
    subplot(1,2,2)    
    hold on
    tt = 0 : Ts : Ts*(length(X_test{i})-1);
    plot(tt, traj_bounds(max_err_ind,:), '--r', 'lineWidth',2)
    for i = 1 : length(X_hat)
        tt = 0 : Ts : Ts*(length(X_test{i})-1);
        plot(tt, vecnorm(X_test{i}-X_hat{i}',2,2))
    end
    xlabel('Time (sec)');
    ylabel('Prediction error $||x-\hat{x}||$');
    title('Test error')
    legend('Maximum error bound')
    
    saveas(fig,'figures/test_fit.png')
    
end