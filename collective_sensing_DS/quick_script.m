%quick_script.m

numRuns = 1;
T = 6000;
mean_err = zeros(numRuns,T);
my_err = mean_err;
wsigma = 0.00;
N = 50;
psi = 0;


for i = 1:numRuns
    i
    [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss(0.1, 'T', T,'N', N, 'psi', psi,'w_sigma',w_sigma);
    mean_err(i,:) = abs(mean(y,2)-psi);
    my_err(i,:) = abs(y(:,3)-psi);
end