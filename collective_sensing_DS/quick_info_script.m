% okay let us try to get out some information theory!!  this is a variation
% on the quick script.


N = 40;
%psi       = 0;
k = 0.1;
verbose   = true;         
timed     = false;         
h         = 1;         
T         = 500;         
sigma     = 0.2;         
wsigma    = 0.05;
tau       = 1;
toplot    = 0;

nbins = 500;
n_cond = 10000;


% to start we need a prior
% normal distribution with mean 0 and std. 1

% let us sample this prior many many times
nsamp = 10000;
psi = normrnd(0,1,[1 nsamp]);

y0s = [];
mean_y0s = zeros(1,nsamp);

for i = 1:nsamp
    y0 = normrnd(psi(i), sigma, [1,N]);
    y0s = [y0s y0];
    mean_y0s(i) = mean(y0);
end


edges = linspace(-5,5,nbins+1);
counts = histc(psi, edges);
counts = counts(1:nbins);
bin_width = 10/nbins;
new_edges = linspace(-5-10*bin_width, 5+10*bin_width, nbins+20+1);
midpoints = edges(1:nbins)+bin_width;
new_midpoints = new_edges(1:nbins+20)+bin_width;
prob = counts/sum(counts);

keep_id = unique(floor(logspace(log10(1),log10(T), 50)));
% 
% ys = cell(1,50);
% mean_ys = ys;

% for j = 51:100
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j} = my_y;
%     mean_ys{j} = mean_y;
% end
%        
%         
% save('bins_51_100_info.mat','ys','mean_ys')
% 
% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 101:150
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j} = my_y;
%     mean_ys{j} = mean_y;
% end
%        
%         
% save('bins_101_150_info.mat','ys','mean_ys')
% 

% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 151:200
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j-150} = my_y(:,keep_id);
%     mean_ys{j-150} = mean_y(:,keep_id);
% end
%        
%         
% save('bins_151_200_info.mat','ys','mean_ys')
% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 201:250
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j-200} = my_y(:,keep_id);
%     mean_ys{j-200} = mean_y(:,keep_id);
% end
%        
%         
% save('bins_201_250_info.mat','ys','mean_ys')
% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 251:300
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j-250} = my_y(:,keep_id);
%     mean_ys{j-250} = mean_y(:,keep_id);
% end
%        
%         
% save('bins_251_300_info.mat','ys','mean_ys')
% 
% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 451:500
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j-450} = my_y(:,keep_id);
%     mean_ys{j-450} = mean_y(:,keep_id);
% end
%        
%         
% save('bins_451_500_info.mat','ys','mean_ys')
% 
% ys = cell(1,50);
% mean_ys = ys;
% 
% for j = 351:400
%     j
%     mean_y = zeros(n_cond, T);
%     my_y = mean_y;
%     psi = midpoints(j);
%     for i = 1:n_cond
%         [~, ~, ~, y] = collective_sensing_DS_loss(k, 'T', T, 'N', N, 'sigma', sigma, 'wsigma', wsigma, 'psi', psi);
%          mean_y(i,:) = mean(y,2);
%          my_y(i,:) = y(:,3);
%     end
%     ys{j-350} = my_y(:,keep_id);
%     mean_ys{j-350} = mean_y(:,keep_id);
% end
%        
%         
% save('bins_351_400_info.mat','ys','mean_ys')
% 
% 
% 
% % %%
% % numRuns = 500;
% % T = 1000;
% % mean_err = zeros(numRuns,T);
% % my_err = mean_err;
% % 
% % 
% % for i = 1:numRuns
% %     i
% %     [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss(0.05, 'T', T);
% %     mean_err(i,:) = abs(mean(y,2));
% %     my_err(i,:) = abs(y(:,3));
% % end