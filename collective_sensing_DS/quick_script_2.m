%quick_script.m
clear y3SAVEneg
clear y3SAVEpos
clear ybarSAVEneg
clear ybarSAVEpos

numRuns = 1000;
T = 100;
mean_err = zeros(numRuns,T);
my_err = mean_err;
w_sigma = 0.05;
sigma = 0.1;
N = 25;
k = 0.1;


%do a bunch of runs with positive psi
psi = .025;

for i = 1:numRuns
    i
    [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss(k, 'T', T,'sigma', sigma, 'w_sigma',w_sigma, 'N', N, 'psi', psi,'toplot',0);
    %mean_err_pos(i,:) = abs(mean(y,2)-psi);
    %my_err_pos(i,:) = abs(y(:,3)-psi);
    y3SAVEpos(i,:) = y(:,3);
    ybarSAVEpos(i,:) = mean(y,2);
end


%do a bunch of runs with negative psi
psi = -.025;

for i = 1:numRuns
    i
    [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss(k, 'T', T,'sigma', sigma,'w_sigma',w_sigma, 'N', N, 'psi', psi,'toplot',0);
    %mean_err_neg(i,:) = abs(mean(y,2)-psi);
    %my_err_neg(i,:) = abs(y(:,3)-psi);
    y3SAVEneg(i,:) = y(:,3);
    ybarSAVEneg(i,:) = mean(y,2);
end

%now bin and compute mutual informations
binBarPos=histc(ybarSAVEpos,-.5:0.01:.5)/numRuns;
binBarNeg=histc(ybarSAVEneg,-.5:0.01:.5)/numRuns;
bin3Pos=histc(y3SAVEpos,-.5:0.01:.5)/numRuns;
bin3Neg=histc(y3SAVEneg,-.5:0.01:.5)/numRuns;

for tct=1:T
    info3(tct) = -sum(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct)).*log2(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct))+eps)) - ...
        ( -0.5*sum(bin3Pos(:,tct).*log2(bin3Pos(:,tct)+eps)) + ...
        -0.5*sum(bin3Neg(:,tct).*log2(bin3Neg(:,tct)+eps)) );
    
    infoBar(tct) = -sum(0.5*(binBarPos(:,tct)+binBarNeg(:,tct)).*log2(0.5*(binBarPos(:,tct)+binBarNeg(:,tct))+eps)) - ...
        ( -0.5*sum(binBarPos(:,tct).*log2(binBarPos(:,tct)+eps)) + ...
        -0.5*sum(binBarNeg(:,tct).*log2(binBarNeg(:,tct)+eps)) );
    
    
    
end

figure()
semilogx(info3)
hold on
semilogx(infoBar,'r')
xlabel('time')
ylabel('MI (bits)')
title('\sigma = 0.1  \sigma_{w} = 0.05  k = 0.1  N = 25')
