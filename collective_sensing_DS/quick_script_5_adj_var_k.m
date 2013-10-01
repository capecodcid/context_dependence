%quick_script.m  for a single parameter scan

addpath(genpath('/Users/andrewhartnett/Documents/MATLAB/context_dependence/'))

clear y3SAVEneg
clear y3SAVEpos
clear ybarSAVEneg
clear ybarSAVEpos
clear info3
clear infoBar

degree = 3;

numRuns = 10000;
T = 100;
mean_err = zeros(numRuns,T);
my_err = mean_err;
w_sigma = 0.05;
N = 100;
sigma= 0.1;
k = 0.1;
%k1 = 0.1;
sig = 0.025;

% % make the graph
A = createRandRegGraph(N, degree);
A = full(A);

% A = ones(N);

%Ns = unique(floor(logspace(log10(1),log10(100),10)));

k1s = linspace(0.05,0.5,10);
%k1s = 0.1;
kct_max = length(k1s);


for kct=1:kct_max
    %N = Ns(kct);
    k1 = k1s(kct);
    
    kct/kct_max
    
    %do a bunch of runs with positive psi
    psi = sig;
    
    for i = 1:numRuns
        [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss_adj_k1(k,k1, A, 'T', T,'w_sigma',w_sigma, 'sigma',sigma,'N', N, 'psi', psi,'toplot',0);
        %mean_err_pos(i,:) = abs(mean(y,2)-psi);
        %my_err_pos(i,:) = abs(y(:,3)-psi);
        y3SAVEpos(i,:) = y(:,1);
        ybarSAVEpos(i,:) = mean(y,2);
    end
    
    
    %do a bunch of runs with negative psi
    psi = -sig;
    
    for i = 1:numRuns
        [best_estimate, best_error, rmserror, y] = collective_sensing_DS_loss_adj_k1(k,k1, A, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'N', N, 'psi', psi,'toplot',0);
        %mean_err_neg(i,:) = abs(mean(y,2)-psi);
        %my_err_neg(i,:) = abs(y(:,3)-psi);
        y3SAVEneg(i,:) = y(:,1);
        ybarSAVEneg(i,:) = mean(y,2);
    end
    
    %now bin and compute mutual informations
    binBarPos=histc(ybarSAVEpos,-.5:0.01:.5)/numRuns;
    binBarNeg=histc(ybarSAVEneg,-.5:0.01:.5)/numRuns;
    bin3Pos=histc(y3SAVEpos,-.5:0.01:.5)/numRuns;
    bin3Neg=histc(y3SAVEneg,-.5:0.01:.5)/numRuns;
    
    for tct=1:T
        info3(kct,tct) = -sum(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct)).*log2(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct))+eps)) - ...
            ( -0.5*sum(bin3Pos(:,tct).*log2(bin3Pos(:,tct)+eps)) + ...
            -0.5*sum(bin3Neg(:,tct).*log2(bin3Neg(:,tct)+eps)) );
        
        infoBar(kct,tct) = -sum(0.5*(binBarPos(:,tct)+binBarNeg(:,tct)).*log2(0.5*(binBarPos(:,tct)+binBarNeg(:,tct))+eps)) - ...
            ( -0.5*sum(binBarPos(:,tct).*log2(binBarPos(:,tct)+eps)) + ...
            -0.5*sum(binBarNeg(:,tct).*log2(binBarNeg(:,tct)+eps)) );
    end
%     figure(1)
%     clf
%     plot(Ns(kct), info3')
%     figure(2)
%     clf
%     plot(Ns(kct), infoBar')
%     drawnow
end


Ns = k1s;
figure()
plot(info3')
hold on
plot(infoBar','r')

figure()
semilogx(info3')
hold on
semilogx(infoBar')


Ts = zeros(1,kct_max);
MImax = Ts;
for i = 1:kct_max
MImax(i) = max(info3(i,:));
Ts(i) = find(info3(i,:)==MImax(i),1);
end


figure()
plot(Ns,infoBar(:,1),'r-*')
hold on
plot(Ns,MImax,'-*')

figure()
plot(Ns,Ts,'r-*')


figure()
plot(MImax,Ts)