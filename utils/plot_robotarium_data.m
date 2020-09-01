clc; clf; clear; close all; addpath('../../koopman_cbf_robotarium_submission/data')

obstacle_avoidance_f = 'obstacle_avoidance_1.mat';
collision_avoidance_f = 'collision_avoidance_6.mat';
collision_obstacle_avoidance_f = 'collision_obstacle_avoidance_12.mat';
f_num = 1;
lw = 2;
ss = 40;
fs = 14;
plot_margin_x = 0.7;
plot_margin_y = 0;

% Plot obstacle avoidance data:
if exist('obstacle_avoidance_f','var') && ~isempty(obstacle_avoidance_f)
    load(obstacle_avoidance_f);
    X = x_data{1};
    backup_val = backup_data{1};
    backup_active = find(backup_val > 1e-3);
    X_init = x_init_data{1};
    X_final = x_final_data{1};
    n_data_pts = size(X,2);
    
    fig = figure(f_num);
    hold on;    
    colors = get(gca,'colororder');
    
    plot(X(1,:),X(2,:),'Color', colors(1,:),'Linewidth',lw)
    plot([X_init(1,1) X_final(1,1)],[X_init(2,1) X_final(2,1)], ':r','Linewidth',lw/2)
    scatter(X(1,backup_active), X(2,backup_active),ss,colors(1,:),'*')
    draw_circle(obs(1),obs(2),r_obs-r_margin/2, 'k')
    
    %axis equal
    title('Koopman CBF single robot obstacle avoidance')
    xlabel('x-coordinate (m)')
    ylabel('y-coordinate (m)')
    legend('Robot trace', 'Desired trajectory', 'Backup controller active','Location','northwest')   
    set(gca,'FontSize',fs)
    
    fig_name = strcat(obstacle_avoidance_f(6:end-4),'.png');
    saveas(fig,strcat('../figures/',fig_name));
    f_num = f_num + 1;
end

% Plot collision avoidance data:
if exist('collision_avoidance_f','var') && ~isempty(collision_avoidance_f)
    load(collision_avoidance_f);
    
    fig = figure(f_num);
    hold on; 
    colors = get(gca,'colororder');
    axis_min = min(min(x_data{1}(1,:)),min(x_data{1}(2,:)));
    axis_max = max(max(x_data{1}(1,:)),max(x_data{1}(2,:)));
    
    legend_plots = [];
    legend_entries = [];
    for i = 1 : length(x_data)
        X = x_data{i};
        backup_val = backup_data{i};
        backup_active = find(backup_val > 2e-1);
        X_init = x_init_data{i};
        X_final = x_final_data{i};
        n_data_pts = size(X,2);

        legend_plots = [legend_plots plot(X(1,:),X(2,:),'Color', colors(i,:),'Linewidth',lw)];
        legend_entries = [legend_entries strcat("Robot ",num2str(i))];
        %legend_plots = [legend_plots scatter(X(1,backup_active), X(2,backup_active),5,colors(i,:),'*')];
        %legend_entries = [legend_entries "CBF active"];
        if i == length(x_data)
            legend_plots = [legend_plots plot([X_init(1,1) X_final(1,1)],[X_init(2,1) X_final(2,1)], ':r','Linewidth',lw/2)];
            legend_entries = [legend_entries "Desired trajectory"];
        else
            plot([X_init(1,1) X_final(1,1)],[X_init(2,1) X_final(2,1)], ':r','Linewidth',lw/2)
        end
        
        axis_min_tmp = min(min(x_data{i}(1,:)),min(x_data{i}(2,:)));
        axis_max_tmp = max(max(x_data{i}(1,:)),max(x_data{i}(2,:)));
        axis_min = min(axis_min, axis_min_tmp);
        axis_max = max(axis_max, axis_max_tmp);
    end
    
    %axis equal
    title('Koopman CBF multi-robot collision avoidance')
    xlabel('x-coordinate (m)')
    ylabel('y-coordinate (m)')
    xlim([axis_min-plot_margin_x axis_max]);
    ylim([axis_min axis_max+plot_margin_y]);
    %reorder_inds = [1 : 2 : 2*length(x_data) 2*length(x_data)+1 2 : 2 : 2*length(x_data)];
    %legend(legend_plots(reorder_inds), legend_entries(reorder_inds), 'Location', 'northwest', 'NumColumns', 2); 
    legend(legend_plots, legend_entries, 'Location', 'northwest', 'NumColumns', 1); 
    set(gca,'FontSize',fs)
    
    fig_name = strcat(collision_avoidance_f(6:end-4),'.png');
    saveas(fig,strcat('../figures/',fig_name));
    f_num = f_num + 1;
end

% Plot collision and obstacle avoidance data:
if exist('collision_obstacle_avoidance_f','var') && ~isempty(collision_obstacle_avoidance_f)
    load(collision_obstacle_avoidance_f);
    
    fig = figure(f_num);
    hold on;    
    colors = get(gca,'colororder');
    axis_min = min(min(x_data{1}(1,:)),min(x_data{1}(2,:)));
    axis_max = max(max(x_data{1}(1,:)),max(x_data{1}(2,:)));
    
    legend_plots = [];
    legend_entries = [];
    n_stop = [1155 1232 1200 1170 1170];
    for i = 1 : length(x_data)
        X = x_data{i};
        backup_val = backup_data{i};
        backup_active = find(backup_val > 2e-1);
        X_init = x_init_data{i};
        X_final = x_final_data{i};
        n_data_pts = size(X,2);

        legend_plots = [legend_plots plot(X(1,1:n_stop(i)),X(2,1:n_stop(i)),'Color', colors(i,:),'Linewidth',lw)];
        legend_entries = [legend_entries strcat("Robot ",num2str(i))];
        %legend_plots = [legend_plots scatter(X(1,backup_active), X(2,backup_active),5,colors(i,:),'*')];
        %legend_entries = [legend_entries "CBF active"];
        if i == length(x_data)
            legend_plots = [legend_plots plot([X_init(1,1) X_final(1,1)],[X_init(2,1) X_final(2,1)], ':r','Linewidth',lw/2)];
            legend_entries = [legend_entries "Desired trajectory"];
        else
            plot([X_init(1,1) X_final(1,1)],[X_init(2,1) X_final(2,1)], ':r','Linewidth',lw/2)
        end
        
        axis_min_tmp = min(min(x_data{i}(1,:)),min(x_data{i}(2,:)));
        axis_max_tmp = max(max(x_data{i}(1,:)),max(x_data{i}(2,:)));
        axis_min = min(axis_min, axis_min_tmp);
        axis_max = max(axis_max, axis_max_tmp);
    end
    draw_circle(obs(1),obs(2),r_obs-r_margin/2, 'k')
    
    %axis equal
    %title('Koopman CBF multi-robot collision and obstacle avoidance')
    xlabel('x-coordinate (m)')
    ylabel('y-coordinate (m)')
    xlim([axis_min-plot_margin_x axis_max]);
    ylim([axis_min axis_max+plot_margin_y]);
    %reorder_inds = [1 : 2 : 2*length(x_data) 2*length(x_data)+1 2 : 2 : 2*length(x_data)];
    %legend(legend_plots(reorder_inds), legend_entries(reorder_inds), 'Location', 'northwest', 'NumColumns', 2);   
    legend(legend_plots, legend_entries, 'Location', 'northwest', 'NumColumns', 1); 
    set(gca,'FontSize',fs)
    
    fig_name = strcat(collision_obstacle_avoidance_f(6:end-4),'.png');
    saveas(fig,strcat('../figures/',fig_name));
    f_num = f_num + 1;
end
        