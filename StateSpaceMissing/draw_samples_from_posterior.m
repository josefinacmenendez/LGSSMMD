function x_l = draw_samples_from_posterior(a, x_k_given_K, sig_sq_i_given_K,N)
L = length(x_k_given_K);
state_space_covariance_vector = a .* sig_sq_i_given_K(2:L);

x_l = zeros(N,L);
for n = 1:N
x_l(n,1) = normrnd( x_k_given_K(1), sqrt( sig_sq_i_given_K(1) ) ) ;
for l = 1:L-1
    x_i_given_K               = x_k_given_K(l);
    x_i_plus_one_given_K      = x_k_given_K(l+1);
    sig_sq_i_plus_one_given_K = sig_sq_i_given_K(l+1);
    sig_sq_i_i_plus_one_given_K = state_space_covariance_vector(l);
    mu_l = x_i_plus_one_given_K + (( sig_sq_i_i_plus_one_given_K / sig_sq_i_plus_one_given_K ) * (x_l(n,l) - x_i_given_K) );
    var_l= sig_sq_i_plus_one_given_K - ( sig_sq_i_i_plus_one_given_K^2 / sig_sq_i_plus_one_given_K );
    x_l(n,l+1) = normrnd( mu_l , sqrt(var_l) ) ;
end
end