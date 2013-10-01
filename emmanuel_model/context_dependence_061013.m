%% let us first build our graph - we are going to use random regular graphs to start
N = 10;
deg = 4;
T = 250000;
y = zeros(T,N);
k1 = 0.0;
k2 = 0.1;
k3 = 1 - k1 - k2;

A = createRandRegGraph(N, deg);
A = full(A)/deg;

XY = rand(N,2);
figure()
gplot(A,XY, '*-')

%% parameters

psi = 0;        % true value - set to zero without loss of generality
sigma = 0.2;       % sampling variance

wsigma = 0.1;    % signaling variance


%% get t=0 samples

% i don't know how to implement the drawing y0 values with a given
% covariance as a function of graph distance.

y0 = normrnd(psi, sigma, [1,N]);
y(1,:) = y0;
best_estimate = mean(y(1,:))

for t = 2:T
    w = normrnd(0, wsigma, [N,N]);
    y(t,:) = k1*y(t-1,:)+(k2)*sum(A.*(repmat(y(t-1,:)',1,N)+ w),1) + k3*y0;
end


figure()
semilogx(sqrt(mean((y.*y)')))
hold on

%% do a second k


y1 = zeros(T,N);
k1 = 0.45;
k2 = 0.1;
k3 = 1-k1-k2;


y1(1,:) = y0;
best_estimate = mean(y1(1,:))

for t = 2:T
    w = normrnd(0, wsigma, [N,N]);
    y1(t,:) = k1*y1(t-1,:)+k2*sum(A.*(repmat(y1(t-1,:)',1,N)+ w),1)+k3*y0;
end

semilogx(sqrt(mean((y1.*y1)')),'r')


%% 
% I should be able to plot information loss with each update.

