function ks = generate_ks(kseed,k_sigma)
    N = length(kseed);
    ks = kseed;
    
    for i = 1:N
        ks(i) = max(normrnd(kseed(i),k_sigma), 0.0);
    end

end