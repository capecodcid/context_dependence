% info_noloss_update.m 
%
% a short function to update the opinions of each agent in our system where
% we don't
%
% Parameters:
%
% N = the number of agents
% T = the total time
% h = the length of a timestep
% psi = the true value
% sigma = the sampling or frozen variance
% w_sigma = the signaling (dynamic) variance
% tau = integration time or 1/k where k is the spring constant
% J = the coupling strength.  Note at J*tau = 1 the system loses stability
% mode = the update rule we want to use
%
%
% Author: Andrew T. Hartnett
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)


function [best_estimate, best_error, rmserror, y, params] = info_noloss_update(varargin)

if (isstruct(varargin)) 
    args= prepareArgs(varargin{1});
else
    args= prepareArgs(varargin);
end

    [   N       ...
        psi     ...
        h       ...
        T       ...
        sigma   ...
        w_sigma  ...
        tau    ...
        J       ...
        mode      ...
        toplot ...
    ] = process_options(args    , ...
    'N'         , 100       , ...
    'psi'       , 0.1          , ...
    'h'         , 1         , ...
    'T'         , 6400         , ...
    'sigma'     , 0.2         , ...
    'w_sigma'    , 0.1         , ...
    'tau'       , 1           , ...
    'J'         , 0.99        , ...
    'mode'      , 'bill' , ...
    'toplot', 0);

    params.psi = psi;
    params.N = N;
    params.h = h;
    params.t = T;
    params.sigma = sigma;
    params.w_sigma = w_sigma;
    params.mode = mode;
    params.J = J;
    params.tau = tau;


mode = 'bill';
N = 10;
h=1;
T=6400;
y = zeros(T,N);


% parameters
psi = 0.1;  % true value - set to zero without loss of generality?
sigma = 0.2;  % sampling variance
w_sigma = 0.1;  % signaling variance
tau=1;  %integration time
J=.999;  %coupling strength. NOTE: at J*tau = 1, the system loses stability


%get t=0 samples
y0 = normrnd(psi, sigma, [1,N]);
y(1,:) = y0;
best_estimate = mean(y(1,:));
best_error=sigma/sqrt(N);

switch mode
    case 'bill'

        for t = 2:T/h
            w = normrnd(0, w_sigma, [1,N]);
            inputs = J/(N-1)*sum(y(t-1,:));
            y(t,:) = y(t-1,:)+(-1/tau*(y(t-1,:)-y0)+inputs-J/(N-1)*y(t-1,:))*h + h^.5*w;
        end
        rmserror=sqrt(mean((y(T/h,:)*(1-J*tau)/J-psi).^2)); 
        
    otherwise
        error('Not a valid mode')
end



if (toplot)
    figure()
    semilogx(mean(y')*(1-J*tau))
    %hold on
    %semilogx([0,T/h+100],[psi/(1-J*tau),psi/(1-J*tau)],'k')
    % i think the point here is that you have to correct for the fact that everything
    % is blowing up as you approach J*tau = 1.

    figure()
    plot(mean(y(1:200,:)')*(1-J*tau)/J)
    % in this figure, you can check how well each individual has approximated
    % the true value
    figure()
    imagesc(y'*(1-J*tau))
    colorbar
end

end
