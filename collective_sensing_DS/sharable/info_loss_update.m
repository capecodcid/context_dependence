% info_loss_update.m
%
% a short function to update the opinions of each agent in our system where
% we bleed out info
%
% Parameters:
%
% ks = vector of learning rates
% psi = true signal
% h = timestep
% T = total time
% sigma = initial sampling sd (quenched noise)
% w_sigma = signalling noise (dynamic noise)
%
%
% Author: Andrew T. Hartnett
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

function [best_estimate, best_error, rmserror, y] = info_loss_update(ks, adj, varargin)

if (isstruct(varargin)) 
    args= prepareArgs(varargin{1});
else
    args= prepareArgs(varargin);
end

    [   psi     ...
        h       ...
        T       ...
        sigma   ...
        w_sigma  ...
        toplot    ...
        mode      ...
        verbose   ...
        meanfield ...
    ] = process_options(args    , ...
    'psi'       , 0          , ...
    'h'         , 1         , ...
    'T'         , 100         , ...
    'sigma'     , 0.2         , ...
    'w_sigma'    , 0.01         , ...
    'toplot'    , 0           , ...
    'mode'      , 'sync' , ...
    'verbose' , 1        , ...
    'meanfield', 0);

    N = length(ks);
    if (size(adj,1) ~= N)
        error('N does not match adj matrix.')
    end
    
    % add error checking to figure out whether or not  N is correct
    
    if (meanfield)
        adj = ones(N);
    end
    
    if (adj == ones(N))  % saves some time
        meanfield = 1;
    end
    
    % t = 0 setup
    y0 = normrnd(psi, sigma, [1,N]);
    y(1,:) = y0;
    best_estimate = mean(y(1,:));
    best_error=sigma/sqrt(N);
    
    
    switch mode
        case 'sync'   % every agent updates every round but in random order
            
            if (meanfield)   % we can speed things up in the meanfield case
                get_inputs = @(y_prev)((ks/N)*sum(y_prev));
            else
                get_inputs = @(y_prev)((ks./sum(adj)).*sum(adj.*repmat(y_prev,N,1)'));
            end
            
            for t = 2:T/h
                
                if (verbose)
                    if mod(t,1000)==0
                        disp(t*h/T)
                    end
                end
                
                w = normrnd(0, w_sigma, [1,N]);
                %inputs = (ks./sum(adj)).*sum(adj.*repmat(y(t-1,:),N,1)');
                inputs = get_inputs(y(t-1,:));
                y(t,:) = (1-ks).*y(t-1,:) + inputs*h + ks.*w;  % do I want to change up the way noise is incorporated
            end
            
            
        case 'async'
            disp('async - not yet coded')

        otherwise
            error('Not a valid mode.')
    end
    
    rmserror=sqrt(mean((y(T/h,:)-psi).^2)); 


end

