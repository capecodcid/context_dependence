%% let us first build our graph - we are going to use random regular graphs to start
N = 10;
deg = 4;
T = 2500;
eta = zeros(T,N);
vis = zeros(T,N);
k = 0.9;

A = createRandRegGraph(N, deg);
A = full(A)/deg;

XY = rand(N,2);
figure()
gplot(A,XY, '*-')

%% in this iteration everyone is required to have a random value on the interval [0,1]

% i am going to start with just a clipped normal distribution although I
% should figure out something more rigorous going forward

psi = 0.5;
sigma = 0.1;

eta(1,:) = normrnd(psi, sigma, [1,N]);
eta(1,eta(1,:)>1) = 1;                      % take care of outliers
eta(1,eta(1,:)<0) = 0;

best_estimate = mean(eta(1,:))

test_vec = rand(1,N);
vis(1,:) = test_vec < eta(1,:);


%% simulate future timesteps

for t = 2:T
    eta(t,:) = k*eta(t-1,:) + (1-k)*(sum(A.*repmat(vis(t-1,:)',1,N),1)-eta(t-1,:));
    test_vec = rand(1,N);
    vis(t,:) = test_vec < eta(t,:);
end

figure()
semilogx(sqrt(mean(((eta-psi).*(eta-psi))')))