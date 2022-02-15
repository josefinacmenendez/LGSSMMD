T = 181;  J_max = 4; sig_e_0 = 2; sig_v_0 = 10; S_0 = 180; t_vec = 1:T;
data_sim = simulate_data(T, J_max, sig_e_0, sig_v_0, S_0);
data_struct = format_data_ss_model(data_sim.Ys, t_vec);
% Estimate initial conditions. 
% If the algorithm doesn't converge, play around with the initial guess for
% the poseterior mode and variance
params=struct('algorithm_params',struct(),'model_params',struct());
params.algorithm_params.('estimate_initial_conditions') = true;
params.algorithm_params.('z_t_given_t_minus_one_0_old') = 220;
params.algorithm_params.('sig_sq_t_given_t_minus_one_0_old') = 500;
params.algorithm_params.('tol')  = 10e-5; 
params.model_params.('sig_sq_e') = 1;
params.model_params.('sig_sq_v') = 1;
model_0 = run_forward_backward_algorithms(data_struct, params);
% Estimate posterior mode and variance given the inidial conditions
params.algorithm_params.('estimate_initial_conditions')= false;
params.algorithm_params.('z_t_given_t_minus_one_0_old') = model_0.z_t_given_t_0;
params.algorithm_params.('sig_sq_t_given_t_minus_0_old')= model_0.sig_sq_t_given_T_0;
params.algorithm_params.tol = 10e-15;
params.model_params.('sig_sq_e') = model_0.sig_sq_e;
params.model_param.s.('sig_sq_v')= model_0.sig_sq_v;
model_1 = run_forward_backward_algorithms(data_struct, params);
% Compute confidence intervals by drawing samples from the posterior distribution
tic
stats_params= struct();
stats_params.('alpha')  = 0.05;    % Confidence level
stats_params.('N')         = 1000; %Number of monte carlo samples
model_stats                   = get_model_stats(model_1,stats_params);
toc
% Plot stats
data_fn = 'testing';
figure_params = struct();
figure_params.('start_idx') = 1;
figure_params.('physiological_data') = {'Heart Rate (BPM)'}; %measured parameter
figure_params.('physiological_str')  = {'heart rate (BPM)'};
figure_params.('phys_id')            = '$Pr(HR_k < HR_j)$';
figure_params.('fig_1_id')           = 'C.';
figure_params.('fig_2_id')           = 'D.';
figure_params.('fn') = strcat(data_fn); %filename to save a figure with the 95% CIs and 
                                                               %a bayesian comparison at each timepoint
get_phys_char_loc_plots(data_struct, model_stats, model_1, figure_params)
