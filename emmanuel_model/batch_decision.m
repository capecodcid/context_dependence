% batch decisions

numruns = 25;
k1s = linspace(0.7,0.89,20);
k3 = 0.1;
k2s = 1-k1s-k3;
deg = 4;
N = 10;
psi = 0.0;
sigma = 1.0;
%wsigmas = linspace(0.1,1,10);
wsigmas = 1.0;

A = createRandRegGraph(N, deg);
A = full(A)/deg;
y0 = normrnd(psi, sigma, [1,N]);

kbs = zeros(length(wsigmas), length(k1s));
kts = kbs;
ks_bs = kbs;
ksprds = kbs;

for w = 1:length(wsigmas)
    for k = 1:length(k1s)
        w
        k
    
        bs = zeros(1,numruns);
        s_bs = bs;
        ts = bs;
        sprds = bs;

        for i = 1:numruns
            [bs(i), s_bs(i), ts(i), sprds(i)] = continuous_iterated_decision(N,deg,4000,k1s(k),k2s(k),k3,sigma,wsigmas(w),0, y0, A);   
        end
    
        kbs(w,k) = mean(bs);
        kts(w,k) = mean(ts);
        ksprds(w,k) = mean(sprds);
        ks_bs(w,k) = mean(s_bs); 
    end
end

%save('deg_4_highw.mat', 'k1s', 'kbs', 'kts', 'ks_bs', 'ksprds');
