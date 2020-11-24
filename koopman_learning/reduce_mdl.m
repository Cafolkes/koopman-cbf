function [K_red, dict_red, C_red] = reduce_mdl(K, func_dict, C)
   threshold = 1e-3;
   
   [~, cols] = find(abs(C)>= threshold);
   used_obs = cols;
   new_obs = cols;
   
   n_obs = 0;
   while n_obs ~= length(used_obs)
       new_obs_tmp = [];
      for i = 1:length(new_obs)
          row = new_obs(i);
          [~, cols] = find(abs(K(row,:)) >= threshold);
          new_obs_tmp = [new_obs_tmp cols(~ismember(cols, used_obs) & ~ismember(cols,new_obs_tmp))];
      end
      new_obs = new_obs_tmp;
      n_obs = length(used_obs);
      used_obs = [used_obs; new_obs_tmp'];
   end
   used_obs = sort(used_obs);
   K_red = K(:,used_obs);
   K_red = K_red(used_obs,:);
   dict_indicator = eye(size(K,1));
   dict_indicator = dict_indicator(used_obs,:);
   dict_red = @(x) dict_indicator*func_dict(x);
   C_red = C(:,used_obs);
end