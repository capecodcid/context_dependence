% collective_sensing_func.m
%
% DS and AH 08-15-2013
%
% parameters:
%
% N = number of agents
% psi = true signal
% h = timestep
% T = total time
% sigma = initial sampling sd (quenched noise)
% w_sigma = signalling noise (dynamic noise)
% tau = -1/tau is the spring constant pulling me back to my original value
% J = coupling strength

function [best_estimate, best_error, rmserror, y] = collective_sensing_DS_evolution_func(ks, adj,varargin)


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
        tau     ...
        toplot    ...
    ] = process_options(args    , ...
    'N'         , 100            , ...
    'psi'       , 0          , ...
    'h'         , 1         , ...
    'T'         , 100         , ...
    'sigma'     , 0.2         , ...
    'w_sigma'    , 0.01         , ...
    'tau'       , 1           , ...
    'toplot'    , 0);

% get t=0 samples
y0 = normrnd(psi, sigma, [1,N]);
y(1,:) = y0;
best_estimate = mean(y(1,:));
best_error=sigma/sqrt(N);

for t = 2:T/h
%     if mod(t,1000)==0
%         time = t*h/T
%     end
    w = normrnd(0, w_sigma, [1,N]);
    %inputs = k/N*sum(y(t-1,:));
    inputs = (ks./sum(adj)).*sum(adj.*repmat(y(t-1,:),N,1)');
    y(t,:) = (1-ks).*y(t-1,:) + inputs*h + ks.*w;  % do I want to change up the way noise is incorporated
end

% for t = 2:T/h
%     w = normrnd(0, w_sigma, [1,N]);
%     inputs = J/(N-1)*sum(y(t-1,:));
%     y(t,:) = y(t-1,:)+(-1/tau*(y(t-1,:)-y0)+inputs-J/(N-1)*y(t-1,:))*h + h^.5*w;
% end

%figure(1)
%clf
%hold on
%plot(mean(y'))
%plot([0,T/h+100],[psi/(1-J*tau),psi/(1-J*tau)],'k')
rmserror=sqrt(mean((y(T/h,:)-psi).^2)); 
% i think the point here is that you have to correct for the fact that everything
% is blowing up as you approach J*tau = 1.

% in this figure, you can check how well each individual has approximated
% the true value
%figure(2)
%imagesc(y'*(1-J*tau))
%colorbar
if (toplot)
    figure()
    plot(abs(mean(y,2)))
    hold on
    plot(abs(y(:,3)),'r')
end

end
