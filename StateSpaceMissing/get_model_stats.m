function model_stats = get_model_stats(model,params)
    model_stats      =   struct();
    N                =   params.N;
    alpha            =   params.alpha;
    A                =   model.A; %state-space smoother;
    sig_sq_t_given_T =   model.sig_sq_t_given_T;
    z_t_given_T      =   model.z_t_given_T;
    uci              =   100*(1-alpha/2);
    lci              =   100*alpha/2;
    T                =   length(z_t_given_T);
    %
    rand_draws = draw_samples_from_posterior(A, (z_t_given_T), (sig_sq_t_given_T),N);
    %
    p_vals_mat = nan(T,T);
    for k = 1:T
        for j = 1:k-1 
            p_vals_mat(k,j) = mean( rand_draws(:,k) < rand_draws(:,j)); 
        end
    end
    %
    model_stats.('p_vals_mat') = p_vals_mat;
    model_stats.('cis')        = prctile(rand_draws,[lci, 50, uci]);
    model_stats.('alpha')    = alpha;
    model_stats.('N')          = N;
    model_stats.('rand_draws') = rand_draws;
end