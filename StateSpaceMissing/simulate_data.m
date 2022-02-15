function data_sim = simulate_data(T,J_max, sig_e_0, sig_v_0, S_0)
Y_sim = zeros(J_max, T);
J_sim = randi(J_max, 1, T);
S_tmp = zeros(1,T-1);
S_tmp(1)= S_0;
for t = 1:T
    S_tmp(1)       = S_0 ;%+ normrnd(0, sig_v_0);
    
    for j = 1:J_max            
        S_tmp(t+1) = S_tmp(t)   + normrnd(0,sig_v_0);
        Y_sim(j,t) = S_tmp(t+1) + normrnd(0, sig_e_0);
    end
end
% choose which ones to keep
for t = 1:T
    J_t = J_sim(t);
    J_nans = J_max - J_t;
    nan_idx= randperm(J_max, J_nans);
    Y_sim(nan_idx,t) = 0;
end
S_sim = sum(Y_sim);
%plot(S_sim)
data_sim = struct('S', S_sim, 'J', J_sim, 'Ys', Y_sim, 'T', length(S_sim));
end
