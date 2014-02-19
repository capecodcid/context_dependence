% fixed_params_clean.m
%
%
% Parameters:
%
% N = # of agents
% signal = true value
% numRuns = # of rounds per generation
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

function [info, infoBar, params] = fixed_param_clean(adj, varargin)
tic
if (isstruct(varargin)) 
    args= prepareArgs(varargin{1});
else
    args= prepareArgs(varargin);
end

    [   signal            ...
        numRuns           ...
        k0_mean           ...
        h                 ...
        T                 ...
        sigma             ...
        w_sigma           ...
        toplot            ...
        mode              ...
        verbose           ...
        computer          ...
        meanfield         ...
    ] = process_options(args          ,     ...
    'signal'      ,       0.025       ,     ...
    'numRuns'     ,       1000        ,     ...
    'k0_mean'     ,       0.2         ,     ...
    'h'           ,       1           ,     ...
    'T'           ,       100         ,     ...
    'sigma'       ,       0.1         ,     ...
    'w_sigma'     ,       0.05        ,     ...
    'toplot'      ,       0           ,     ...
    'mode'        ,       'sync'      ,     ...
    'verbose'     ,       1           ,     ...
    'computer'    ,       'andrew'    ,     ...
    'meanfield'   ,       0           );     % lets you override adj

    params.adj = adj;
    params.signal = signal;
    params.k0_mean = k0_mean;
    params.h = h;
    params.t = T;
    params.sigma = sigma;
    params.w_sigma = w_sigma;
    params.mode = mode;

% 
% signal = 0.035;
% numRuns = 1000;
%     k0_mean     =       0.2;
%     k0_sigma    =       0.05;
%     mutation    =       0.01 ;
%     h          = 1  ;
%     T           =       100   ;
%     sigma       =       0.1   ;
%     w_sigma     =       0.05    ;
%     target_info =       0.65      ;     % !ath
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
    
% preallocate some of these for speed??  !ath
N = size(adj,1);

mean_err = zeros(numRuns,T);
my_err = mean_err;

ks = ones(1,N)*k0_mean;
    
    %do a bunch of runs with positive psi
    psi = signal;
    
    for i = 1:numRuns
        [~, ~, ~, y] = info_loss_update(ks, adj, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'psi', psi,'toplot',0, 'meanfield',meanfield);
        ySAVEpos(i,:,:) = y;
        ybarSAVEpos(i,:) = mean(y,2);
    end
    
    
    %do a bunch of runs with negative psi
    psi = -signal;
    
    for i = 1:numRuns
        [~, ~, ~, y] = info_loss_update(ks, adj, 'T', T,'w_sigma',w_sigma, 'sigma',sigma, 'psi', psi,'toplot',0,'meanfield',meanfield);
        ySAVEneg(i,:,:) = y;
        ybarSAVEneg(i,:) = mean(y,2);
    end
    
    %now bin and compute mutual informations
    binBarPos=histc(ybarSAVEpos,-.5:0.01:.5)/numRuns;
    binBarNeg=histc(ybarSAVEneg,-.5:0.01:.5)/numRuns;
    
    binPos=histc(ySAVEpos,-.5:0.01:.5)/numRuns;
    binNeg=histc(ySAVEneg,-.5:0.01:.5)/numRuns;
   
    
    for tct=1:T
        for id = 1:N
            info(tct, id) = -sum(0.5*(binPos(:,tct, id)+binNeg(:,tct, id)).*log2(0.5*(binPos(:,tct, id)+binNeg(:,tct, id))+eps)) - ...
                ( -0.5*sum(binPos(:,tct, id).*log2(binPos(:,tct,id)+eps)) + ...
                -0.5*sum(binNeg(:,tct, id).*log2(binNeg(:,tct, id)+eps)) );
        end
        
        infoBar(tct) = -sum(0.5*(binBarPos(:,tct)+binBarNeg(:,tct)).*log2(0.5*(binBarPos(:,tct)+binBarNeg(:,tct))+eps)) - ...
            ( -0.5*sum(binBarPos(:,tct).*log2(binBarPos(:,tct)+eps)) + ...
            -0.5*sum(binBarNeg(:,tct).*log2(binBarNeg(:,tct)+eps)) );
    end
    
    info = squeeze(info);
%     Ts = zeros(1,N);
%     MImax = Ts;
%     for i = 1:N
%         MImax(i) = max(info(:,i));
%         Ts(i) = find(info(:,i)==MImax(i),1);
%     end
toc  
end

