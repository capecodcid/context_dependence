%tradeoff_script.m

N = 50;
ks_to_run = [0.05:0.05:0.95];
adj = ones(N);

% preallocate storage


for i = 1:length(ks_to_run)
    k = ks_to_run(i);
    [info, ~, ~] = fixed_param_clean(adj,'k0_mean',k);
    infos(i,:) = mean(info,2);
end

%% let us plot our infos up

% build a color map
time = 1:size(infos,2);

cmap = colormap(lbmap(size(infos,1),'RedBlue'));


figure()
semilogx(time,infos(1,:), 'color', cmap(1,:))
hold on

for i = 2:size(infos,1)
    semilogx(time,infos(i,:),'color', cmap(i,:));
end

xlabel('timesteps')
ylabel('MI bits')

%% next we want to look at for each cutoff who gets there first

cutoffs = linspace(0.25,0.7,15);
Ts = zeros(1, length(cutoffs));
best_k = Ts;

for i = 1:length(cutoffs)
    cutoff = cutoffs(i);
    t_ = size(infos,2);
    k_ = 0;
    for j = 1:size(infos,1)
        tmp = find(infos(j,:)>=cutoff,1);
        if (~isempty(tmp) && tmp < t_) 
            t_ = tmp;
            k_ = ks_to_run(j);
        end
    end 
    Ts(i) = t_;
    best_k(i) = k_;
end

figure()
plot(cutoffs, best_k, 'o')
xlabel('cutoff bits')
ylabel('best k value')

figure()
semilogy(cutoffs, Ts, 'ro')
xlabel('cutoff (bits)')
ylabel('shortest time')
