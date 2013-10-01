%quick_script.m  for a single parameter scan

addpath(genpath('/Users/andrewhartnett/Documents/MATLAB/context_dependence/'))

clear y3SAVEneg
clear y3SAVEpos
clear ybarSAVEneg
clear ybarSAVEpos
clear info3
clear infoBar

clear info
clear ySAVEneg
clear ySAVEpos

degree = 3;

numRuns = 1000;     % per generation per outcome
T = 150;                % per discussion
mean_err = zeros(numRuns,T);
my_err = mean_err;
w_sigma = 0.05;    % dynamic (talking) noise
N = 100;
sigma= 0.1;         % quenched noise

k_sigma = 0.05;    % sd for variablity in ks
mutation = 0.01;    % mutation rate

signal = 0.035;
k0 = 0.1;  % starting k
k0 = ones(N,1)*k0;

ks = generate_ks(k0, k_sigma);
ks = ks';

generations = 2;

for gen=1:generations
    %N = Ns(kct);
    %k1 = k1s(kct);
    
    gen/generations
    
    %do a bunch of runs with positive psi
    psi = signal;
    
    for i = 1:numRuns
        %A = createRandRegGraph(N, degree);  % generate a new graph each run
        %A = full(A);
        A = ones(N);
        [~, ~, ~, y] = collective_sensing_DS_evolution_func_async(ks, A, 'T', T,'w_sigma',w_sigma, 'sigma',sigma,'N', N, 'psi', psi,'toplot',0);
        %mean_err_pos(i,:) = abs(mean(y,2)-psi);
        %my_err_pos(i,:) = abs(y(:,3)-psi);
        %y3SAVEpos(i,:) = y(:,1);
        ySAVEpos(i,:,:) = y;
        ybarSAVEpos(i,:) = mean(y,2);
    end
    
    
    %do a bunch of runs with negative psi
    psi = -signal;
    
    for i = 1:numRuns
        %A = createRandRegGraph(N, degree);  % generate a new graph each run
        %A = full(A);
        A = ones(N);
        [~, ~, ~, y] = collective_sensing_DS_evolution_func_async(ks, A, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'N', N, 'psi', psi,'toplot',0);
        %mean_err_neg(i,:) = abs(mean(y,2)-psi);
        %my_err_neg(i,:) = abs(y(:,3)-psi);
        %y3SAVEneg(i,:) = y(:,1);
        ySAVEneg(i,:,:) = y;
        ybarSAVEneg(i,:) = mean(y,2);
    end
    
    %now bin and compute mutual informations
    binBarPos=histc(ybarSAVEpos,-.5:0.01:.5)/numRuns;
    binBarNeg=histc(ybarSAVEneg,-.5:0.01:.5)/numRuns;
    
    binPos=histc(ySAVEpos,-.5:0.01:.5)/numRuns;
    binNeg=histc(ySAVEneg,-.5:0.01:.5)/numRuns;
   
    %bin3Pos=histc(y3SAVEpos,-.5:0.01:.5)/numRuns;
    %bin3Neg=histc(y3SAVEneg,-.5:0.01:.5)/numRuns;
    
    for tct=1:T
        for id = 1:N
            info(tct, id) = -sum(0.5*(binPos(:,tct, id)+binNeg(:,tct, id)).*log2(0.5*(binPos(:,tct, id)+binNeg(:,tct, id))+eps)) - ...
                ( -0.5*sum(binPos(:,tct, id).*log2(binPos(:,tct,id)+eps)) + ...
                -0.5*sum(binNeg(:,tct, id).*log2(binNeg(:,tct, id)+eps)) );
        end
        
%         info3(gen,tct) = -sum(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct)).*log2(0.5*(bin3Pos(:,tct)+bin3Neg(:,tct))+eps)) - ...
%             ( -0.5*sum(bin3Pos(:,tct).*log2(bin3Pos(:,tct)+eps)) + ...
%             -0.5*sum(bin3Neg(:,tct).*log2(bin3Neg(:,tct)+eps)) );
        
        infoBar(gen,tct) = -sum(0.5*(binBarPos(:,tct)+binBarNeg(:,tct)).*log2(0.5*(binBarPos(:,tct)+binBarNeg(:,tct))+eps)) - ...
            ( -0.5*sum(binBarPos(:,tct).*log2(binBarPos(:,tct)+eps)) + ...
            -0.5*sum(binBarNeg(:,tct).*log2(binBarNeg(:,tct)+eps)) );
    end
    info = squeeze(info);
    
    k_record(gen,:) = ks; 
    
    fitness = sum(info);
    
    h = figure('visible','off')
    plot(ks,fitness,'o')
    xlim([0,0.35])
    ylim([0,170])
    saveas(h,strcat('async_generation_meanfield_small_mute_',num2str(gen),'.png'))
    
    
    
    ks = randsample(ks,100,true,fitness);
    
    
    mutate_id = find(rand(1,100) < mutation);
    for j = 1:length(mutate_id)
        ks(mutate_id(j)) = max(normrnd(ks(mutate_id(j)),k_sigma),0.0);
    end
    
   
end

%info = squeeze(info);
% check above - fix below

% 
% figure()
% plot(info3')
% %hold on
% %plot(infoBar','r')
% 
% figure()
% semilogx(info3')
% %hold on
% %semilogx(infoBar')
% 
% 
% Ts = zeros(1,kct_max);
% MImax = Ts;
% for i = 1:kct_max
% MImax(i) = max(info3(i,:));
% Ts(i) = find(info3(i,:)==MImax(i),1);
% end
% 
% 
% figure()
% plot(Ns,infoBar(:,1),'r-*')
% hold on
% plot(Ns,MImax,'-*')
% 
% figure()
% plot(Ns,Ts,'r-*')
% 
% 
% figure()
% plot(MImax,Ts)