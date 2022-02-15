function data_struct = format_data_ss_model( data, t_vec)
% This function takes as input a N by T matrix with data that may contain
% missing observations and returns a struct with the sum of the
% observations, the length of the timeseries, and the number of
% observations at each timepoint.
% Inputs: data: <double> : A N by T matrix with data points
%             t_vec:<double> : A 1 by T vector with times at which data was
%             collected
% Output: 
%               data_struct: A struct with the following fields
%                                   S:  <double>:   A 1 by T vector with
%                                   the sum of observations across trials
%                                    Js: <double>:  A 1 by T vector with
%                                    the number of observations at each
%                                    timepoint
%                                   Ys: <double>: A N by T matrix with all
%                                   observations
%                                   T: <integer>: The length of the
%                                   timeseries
    [r_idxs,~] = find(~isnan(data));  
    trials_with_data  = unique(r_idxs);
    Ys        = data(trials_with_data,:); %Extract trials with data
    Ys(isnan(Ys)) = 0; %Set nan's to 0 
    Js        = sum( ~isnan(data) ,1); % extract row (r) indexes and column (c) indexes with missing data (nan)
    S         = nansum(data);             % compute sum of all values of data for J sessions
    T         = length(S);                     % compute length of the time-series (T time points)
    data_struct = struct('T', T, 'S', S, 'Ys', Ys, 'J', Js, 'time_vec', t_vec);
end