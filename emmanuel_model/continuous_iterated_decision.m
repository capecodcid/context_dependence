function [best, scale_best, time, spread] = continuous_iterated_decision(N, deg, tmax, k1, k2, k3, sigma, wsigma, fig, y0, A)
    
    psi = 0; %w/o loss of generality
% set up
%     A = createRandRegGraph(N, deg);
%     A = full(A)/deg;
    
    k = k1+k2+k3;
    k1 = k1/k;
    k2 = k2/k;
    k3 = k3/k;
    
    XY = rand(N,2);
    %figure()
    %gplot(A,XY,'*-')
    
    y = zeros(tmax,N);
    
%     y0 = normrnd(psi, sigma, [1,N]);
    best_estimate = mean(y0);
    
    
    y(1,:) = y0;
    for t = 2:tmax
        w = normrnd(0, wsigma, [N,N]);
        y(t,:) = k1*y(t-1,:)+(k2)*sum(A.*(repmat(y(t-1,:)',1,N)+w),1) + k3*y0;
    end
    
    if (fig)
    
        figure()
        hold on
        for i = 1:N
            semilogx(y(:,i),'b-')
        end
        %hline(best_estimate,'k-')
        %hline(-best_estimate,'k-')
        set(gca,'xscale','log')
        xlabel('timesteps','FontSize',14);
        ylabel('estimate','FontSize',14);
        hold off
    
    
        figure()
        loglog(smooth(var(y'),10))
    end
    
    th = 2*mean(var(y(end-500:end,:)'));
    time = find(var(y')<th,1);
    best = sqrt(mean(y(time,:).^2));
    scale_best = best/sqrt(best_estimate^2);
    spread = th/2;
    
end
