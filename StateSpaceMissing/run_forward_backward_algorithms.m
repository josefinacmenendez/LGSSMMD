function model = run_forward_backward_algorithms(data, params)
% This function computes the posterior mode and variance for a linear gaussian state-space model
% Inputs: data: a struct with the following fields:
%		S:					<double>: A 1 by T vector with the sum of observations (across time)
%		Ys:					<double>: A J by T matrix with all observations
%								(J can be number of sessions/subjects, T is the length of the timeseries)
%		J:					<integer>: a 1 by T vector with the number of observations at each timepoint
%	  params: a struct with the following fields:
%		estimate_initial_conditions: 		<bool>:   if true, flips the data, and estimates the initial guess 
%				     				  for the posterior mode and variance
%		z_t_given_t_minus_one_0_old: 		<double>: initial guess for the posterior mode
%		sig_sq_t_given_t_minus_one_0_old:	<double>: initial guess for the posterior variance
%		tol:					<double>: tolerance (for convergence)
%		sig_sq_e:				<double>: initial guess for the observation variance
%		sig_sq_v:				<double>: initial guess for the state process variance

check_convergence = @(prev_param_estimate, curr_param_estimate, tol) abs((prev_param_estimate-curr_param_estimate)/prev_param_estimate)<tol;
log_likelihood    = @(sig_sq_e, sig_sq_v, Ys, z_t_given_T, T,J,W_t, W_t_t_minus_one, W_t_minus_one)...
                        (-0.5*( (2*pi*sig_sq_e*sum(J)) + (1/sig_sq_e) * ...
                        sum( sum(Ys(:,2:end).^2) + (-2* z_t_given_T(2:end)...
                        .* sum(Ys(:,2:end))) + (J(2:end).*W_t)) ...
                        +(T-1)*log(2*pi*sig_sq_v) + (1/sig_sq_v)*...
                        sum( W_t - (2 * W_t_t_minus_one ) + W_t_minus_one )));
compute_sig_sq_e    = @(W_t, z_t_given_T, S, J, Ys) (1/(sum(J))).* sum( sum(Ys(:,2:end).^2) + (-2* z_t_given_T(2:end) .* sum(Ys(:,2:end))) + (J(2:end).*W_t) );
compute_sig_sq_v    = @(W_t, W_t_t_minus_one, W_t_minus_one, T, J) ...
                    (1/(T-1))    * sum( W_t - (2 * W_t_t_minus_one ) + W_t_minus_one );
%
S = data.S;
Ys= data.Ys;
J = data.J;
T = data.T;
estimate_initial_conditions = params.algorithm_params.estimate_initial_conditions;
z_t_given_t_minus_one_0_old = params.algorithm_params.z_t_given_t_minus_one_0_old;
sig_sq_t_given_t_minus_one_0_old = params.algorithm_params.z_t_given_t_minus_one_0_old;
tol = params.algorithm_params.tol;
%
sig_sq_e = params.model_params.sig_sq_e;
sig_sq_v = params.model_params.sig_sq_v;
%
if estimate_initial_conditions
    S = flip(S);
    Ys= flip(Ys,2);
    J = flip(J);
end
%
converged = false;
% initialize parameters
n               = 1;
dev_old         = 1;
sig_sq_e_old    = 1;
sig_sq_v_old    = 1;
%
lls             = {};
%
while ~converged
    sig_sq_t_given_t_minus_one    = zeros(1,T);
    z_t_given_t_minus_one         = zeros(1,T);
    z_t_given_t                   = zeros(1,T);
    sig_sq_t_given_t              = zeros(1,T);
    z_t_given_t(1)                = z_t_given_t_minus_one_0_old;
    sig_sq_t_given_t(1)           = sig_sq_t_given_t_minus_one_0_old;
    % OK
    for t = 2:T
        % one-step prediction
        z_t_given_t_minus_one(t)      = z_t_given_t(t-1);                   %one-step prediction mode;
        sig_sq_t_given_t_minus_one(t) = sig_sq_t_given_t(t-1) + sig_sq_v;   %one-step prediction variance
        % posterior mode and variance
        z_t_given_t(t)                = ((z_t_given_t_minus_one(t) * sig_sq_e) + (S(t)*sig_sq_t_given_t_minus_one(t)))/...
                                        (sig_sq_e + (J(t) * sig_sq_t_given_t_minus_one(t)));
        sig_sq_t_given_t(t)           = (sig_sq_e * sig_sq_t_given_t_minus_one(t)) / ( sig_sq_e + ( J(t) * sig_sq_t_given_t_minus_one(t) ) );
    end    
    % FIS
    z_t_given_T             = zeros(1,T);
    sig_sq_t_given_T        = zeros(1,T);
    % initial guess
    z_t_given_T(end)        = z_t_given_t(end);
    sig_sq_t_given_T(end)   = sig_sq_t_given_t(end);
    A                       = zeros(1,T-1);
    for t = T-1: -1 : 1
        A(t)                =   sig_sq_t_given_t(t) / sig_sq_t_given_t_minus_one(t+1);
        z_t_given_T(t)      =   z_t_given_t(t) + ( A(t)*(z_t_given_T(t+1) - z_t_given_t_minus_one(t+1)) );
        sig_sq_t_given_T(t) =   sig_sq_t_given_t(t) + ( A(t)^2* (sig_sq_t_given_T(t+1)-sig_sq_t_given_t_minus_one(t+1)) );
    end   
    % state-space covariance algorithm
    sig_t_minus_one_t_given_T= A .* sig_sq_t_given_T(2:end);
    % Update variance
    S_j                      = sum(J(2:end));
    W_t                      = sig_sq_t_given_T(2:end)   + (z_t_given_T(2:end).^2);
    W_t_t_minus_one          = sig_t_minus_one_t_given_T + (z_t_given_T(1:end-1) .* z_t_given_T(2:end));
    W_t_minus_one            = sig_sq_t_given_T(1:end-1) + (z_t_given_T(1:end-1) .^2);
    % Update variances
    sig_sq_e = compute_sig_sq_e(W_t, z_t_given_T, S, J, Ys);
    sig_sq_v = compute_sig_sq_v(W_t, W_t_t_minus_one, W_t_minus_one, T, J);
    % Check for convergence   
    ll        = log_likelihood(sig_sq_e, sig_sq_v, Ys, z_t_given_T, T,J,W_t, W_t_t_minus_one, W_t_minus_one);
    lls{n}    = ll;
    dev_new   = -2*ll;
    if estimate_initial_conditions
        convergence =   check_convergence(sig_sq_e_old, sig_sq_e, tol) && ...
                        check_convergence(sig_sq_v_old, sig_sq_v, tol) && ...
                        check_convergence(z_t_given_t_minus_one_0_old, z_t_given_T(end), tol) && ...
                        check_convergence(sig_sq_t_given_t_minus_one_0_old, sig_sq_t_given_T(end), tol) && ...
                        check_convergence(dev_old, dev_new, tol);
    else
        convergence =   check_convergence(sig_sq_e_old, sig_sq_e, tol) && ...
                        check_convergence(sig_sq_v_old, sig_sq_v, tol) && ...
                        check_convergence(dev_old, dev_new, tol);
    end
    if convergence
        fprintf('Converged in %s iterations.\n', string(n));
        converged = true;
        if estimate_initial_conditions
            model = struct('z_t_given_t_0',z_t_given_T(end), ...
			   'sig_sq_t_given_T_0',sig_sq_t_given_T(end), ...
			   'sig_sq_e', sig_sq_e, ...
		           'sig_sq_v', sig_sq_v);
        else
            model = struct('z_t_given_T', z_t_given_T, ...
                            'sig_sq_t_given_T', sig_sq_t_given_T, ...
                            'sig_sq_e', sig_sq_e, ...
                            'sig_sq_v', sig_sq_v, ...
                            'A'       , A);
        end
    else
        n = n + 1;
        if estimate_initial_conditions
            z_t_given_t_minus_one_0_old       = z_t_given_T(end);
            sig_sq_t_given_t_minus_one_0_old  = sig_sq_t_given_T(end);
            %
        end
        dev_old         = dev_new;
        sig_sq_e_old    = sig_sq_e;
        sig_sq_v_old    = sig_sq_v;
    end
end
end
