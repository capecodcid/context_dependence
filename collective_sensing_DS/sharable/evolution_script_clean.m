% evolution_script_clean.m
%
% here we allow the parameters ks to evolve based upon some fitness
% function
%
% Parameters:
%
% N = # of agents
% signal = true value
% numRuns = # of rounds per generation
% generations = # of generations
% k0_mean = mean starting k
% k0_sigma = spread of starting ks
% mutation = mutation rate
% h = timestep
% T = time per discussion
% sigma = initial sampling sd (quenched noise)
% w_sigma = signalling noise (dynamic noise)
% target_info = how many bits our agents need to thrive              % !ath
% computer = which computer this is being run on
% meanfield = should we do the meanfield speedup
%
% OUTPUT:
% info gives the mutual information for each individual at each time point
% for the last generation it is T-by-N
%
% infoBar gives the mutual information contained by the mean of all agents
% at each timepoint in each generation .. it is gen-by-T
%
% k_record gives the k of each agent in each generation gen-by-N
%
% mean_k and sd_k give mean and stdev by generation
%
% Author: Andrew T. Hartnett
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

function [info, infoBar, k_record, params] = evolution_script_clean(adj, varargin)

if (isstruct(varargin)) 
    args= prepareArgs(varargin{1});
else
    args= prepareArgs(varargin);
end

    [   signal            ...
        numRuns           ...
        generations       ...
        k0_mean           ...
        k0_sigma          ...
        mutation          ...
        h                 ...
        T                 ...
        sigma             ...
        w_sigma           ...
        target_info       ...          % !ath
        toplot            ...
        mode              ...
        verbose           ...
        computer          ...
        meanfield         ...
    ] = process_options(args          ,     ...
    'signal'      ,       0.035       ,     ...
    'numRuns'     ,       1000        ,     ...
    'generations' ,       20          ,     ...
    'k0_mean'     ,       0.2         ,     ...
    'k0_sigma'    ,       0.05        ,     ...
    'mutation'    ,       0.01        ,     ...
    'h'           ,       1           ,     ...
    'T'           ,       100         ,     ...
    'sigma'       ,       0.1         ,     ...
    'w_sigma'     ,       0.05        ,     ...
    'target_info' ,       0.85        ,     ...       % !ath
    'toplot'      ,       0           ,     ...
    'mode'        ,       'sync'      ,     ...
    'verbose'     ,       1           ,     ...
    'computer'    ,       'andrew'    ,     ...
    'meanfield'   ,       0           );     % lets you override adj

    params.adj = adj;
    params.signal = signal;
    params.k0_mean = k0_mean;
    params.k0_sigma = k0_sigma;
    params.mutation = mutation;
    params.h = h;
    params.t = T;
    params.sigma = sigma;
    params.w_sigma = w_sigma;
    params.target_info = target_info;
    params.mode = mode;


% 
% signal = 0.035;
% numRuns = 1000;
% generations = 10;
%     k0_mean     =       0.2;
%     k0_sigma    =       0.05;
%     mutation    =       0.01 ;
%     h           =       1  ;
%     T           =       100   ;
%     sigma       =       0.1   ;
%     w_sigma     =       0.05    ;
%     target_info =       0.45      ;     % !ath
%     toplot    =      0           ;
%     mode      =        'sync'    ;
%     verbose   =         1        ;
%     computer  =         'andrew'  ;
%     meanfield =         0           ;
    
switch computer
    case 'andrew' 
        addpath(genpath('/Users/andrewhartnett/Documents/MATLAB/context_dependence/'));
        addpath(genpath('/Users/andrewhartnett/Documents/MATLAB/context_dependence/collective_sensing_DS/'))
        addpath(genpath('/Users/andrewhartnett/Documents/MATLAB/context_dependence/randRegGraph/randRegGraph/'))

    case 'hwuimac'
        addpath(genpath('/Users/ahartnet/Documents/git/context_dependence/'));
        addpath(genpath('/Users/ahartnet/Documents/git/context_dependence/collective_sensing_DS/'))
        addpath(genpath('/Users/ahartnet/Documents/git/context_dependence/randRegGraph/randRegGraph/'))
        
    otherwise
        error('not a valid computer choice')
end
    
clear y3SAVEneg
clear y3SAVEpos
clear ybarSAVEneg
clear ybarSAVEpos
clear info3
clear infoBar

clear info
clear ySAVEneg
clear ySAVEpos

N = size(adj,1);

mean_err = zeros(numRuns,T);
my_err = mean_err;

k0_mean = ones(N,1)*k0_mean;

ks = generate_ks(k0_mean, k0_sigma);
ks = ks';

for gen=1:generations
    tic
    gen/generations
    
    %do a bunch of runs with positive psi
    psi = signal;
    
    for i = 1:numRuns
        [~, ~, ~, y] = info_loss_update(ks, adj, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'psi', psi,'toplot',0, 'meanfield',meanfield);
        %mean_err_pos(i,:) = abs(mean(y,2)-psi);
        %my_err_pos(i,:) = abs(y(:,3)-psi);
        %y3SAVEpos(i,:) = y(:,1);
        ySAVEpos(i,:,:) = y;
        ybarSAVEpos(i,:) = mean(y,2);
    end
    
    
    %do a bunch of runs with negative psi
    psi = -signal;
    
    for i = 1:numRuns
        [~, ~, ~, y] = info_loss_update(ks, adj, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'psi', psi,'toplot',0,'meanfield',meanfield);
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
    
    fitness = zeros(N,1);
    for i = 1:N
        tmp = find(info(:,i)>target_info,1);
        if ~isempty(tmp)
            fitness(i) = 1/tmp;
        else
            fitness(i) = 0;
        end
    end
    [donotuse, ix]=sort(fitness);
    [donotuse, yourranks]=sort(ix);
    fitness = yourranks;
    fitness = fitness/sum(fitness);
    
    %fitness = sum(info);load
    %[donotuse, ix]=sort(fitness);
    %[donotuse, yourranks]=sort(ix);
    %fitness = yourranks;
    
    
    Ts = zeros(1,N);
    MImax = Ts;
    for i = 1:N
        MImax(i) = max(info(:,i));
        Ts(i) = find(info(:,i)==MImax(i),1);
    end
    
% %     h = figure('visible','off')
%      figure()
%      plot(ks,fitness,'o')
%      xlim([0,1])
%      ylim([0,0.04])
% %     saveas(h,strcat('async_k0_01_meanfield_fit_05_speedrank_v2_',num2str(gen),'.png'))
% %     
%     
    
    ks = randsample(ks,N,true,fitness);    
    
    mutate_id = find(rand(1,N) < mutation);
    for j = 1:length(mutate_id)
        ks(mutate_id(j)) = max(normrnd(ks(mutate_id(j)),k0_sigma),0.0);
    end   
toc   
end

end
% figure()
% plot(info(:,3)')
% hold on
% plot(infoBar(end,:)','r')
% 
% 
% Ts = zeros(1,N);
% MImax = Ts;
% for i = 1:N
%     MImax(i) = max(info(:,i));
%     Ts(i) = find(info(:,i)==MImax(i),1);
% end
% 
% Ns = [1:N];
% 
% figure()
% plot(Ns,MImax,'-*')
% 
% figure()
% plot(Ns,Ts,'r-*')
%  
%  
% figure()
% plot(MImax,Ts)
