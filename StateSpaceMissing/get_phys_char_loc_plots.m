function get_phys_char_loc_plots(data, model_stats, model, figure_params)
TS = 12;% title font size
LS = 12;% legend font size
CS = 10;% caption font size
XS = 8; %tick font size
set(0, 'defaulttextinterpreter', 'latex')
J             = data.J;            % number of observations per time
S            = data.S;            % sum of observations for each time point
time_vec= data.time_vec;     % time vector
n             = size(data.Ys,1);   % number of sessions used
%
cis         = model_stats.cis;          % confidence intervals
p_vals_mat  = model_stats.p_vals_mat;   % minute-by-minute comparison
z_t_given_T = model.z_t_given_T;        % posterior mode
alpha       = model_stats.alpha;        % confidence level
N           = model_stats.N;            % number of monte carlo samples
%
start_idx   = figure_params.start_idx;  % start plotting at start index
phys_metric = figure_params.physiological_data{:}; % parameter estimated
phys_str    = figure_params.physiological_str{:};  % parameter estimated
phys_id     = figure_params.phys_id;
fig_1_id    = figure_params.fig_1_id;
fig_2_id    = figure_params.fig_2_id;
%
[sig_i, sig_j] = find(vpa(p_vals_mat >= 1-alpha));
min_param_t    = time_vec(min(sig_i));
fig_handle = figure('Renderer','Painters');
fig_handle.PaperUnits = 'inches';
fig_handle.PaperSize = [16 14];
p  = subplot(4,3,[1 2 3 4 5 6]); %[0.05 0.3 0.75]
h  = scatter(time_vec(start_idx:end), S(start_idx:end)./J(start_idx:end),'filled','MarkerFaceColor',[0.05 0.3 0.75],'MarkerFaceAlpha',0.5); hold on;
h1 = plot(time_vec(start_idx:end),(cis(1,start_idx:end)),'-r','linewidth',1.5);hold on;
h2 = plot(time_vec(start_idx:end),(z_t_given_T(start_idx:end)),'-b','linewidth',1.5);hold on;
h3 = plot(time_vec(start_idx:end),(cis(3,start_idx:end)),'r','linewidth',1.5);hold on;
h4 = xline(0 ,'linewidth',1.5,'color',[0.01 0.01 0.01]);
h5 = xline(60,'linewidth',1.5,'color',[0.01 0.01 0.01]);
%
xl = get(gca,'XLabel');
yl = get(gca, 'YLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', XS)
set(xl, 'FontSize', LS);
yAX = get(gca,'YAxis');
set(yAX, 'FontSize', XS);
set(yl, 'FontSize', LS);
%ylim([0 220]);
legend([h,h2, h1, h4], {'Empirical mean','Posterior mean',...
                                    strcat(string(100*(1-alpha)),{' '},...
                                    '% Confidence Intervals'),'Time of propofol infusion'},...
                                    'location','northeast','fontsize', LS);
fig_1_id_t = title(fig_1_id,'fontsize',TS,'fontweight','bold');
fig_1_id_t.Position = [time_vec(start_idx) fig_1_id_t.Position(2:end)];
%
f       = subplot(4,3,[7 10]);
dim     = get(f,'position');
dim(2)  = 0.1;
str =  strcat(fig_1_id,' Blue dots show the mean', {' '} , phys_str, ...
        {' '} ,' for N = ',int2str(n),' sessions. ', ...
        {' '} , 'EM algorithm estimates of the ', {' '} , phys_str, ...
        {' '} ,' (blue line) and its ', {' '} ,string(100*(1-alpha)), ...
        {' '} , '\% confidence intervals (red lines) estimated ', ...
        {' '} , 'using ', {' '} , int2str(N), ...
        {' '} , ' monte carlo samples.  ',...
        {' '} , fig_2_id,{' '},'Minute-by-minute comparison of ',...
        {' '} , 'the probability that the', {' '} , phys_str,...
        {' '} , 'at time k is lower than the ', {' '} , phys_str , ...
        {' '} , ' at time j, ',...
        {' '} , 'where k = 1, $\dots$, 180 and j = 2, $\dots$, 181. Comparisons ',...
        {' '} , 'for which this probability is $>= $', {' '} , string(1-alpha) , ...
        {' '} , ' are shown in red. The earliest ',...
        {' '} , 'statistically significant difference in ', ...
        {' '} , phys_str , {' '} ,' was obtained at ', ...
        {' '} , num2str(min_param_t), {' '} ,' minutes after the propofol infusion.');
annotation('textbox',dim,'String',strjust(str),'interpreter','latex',...
    'fontsize',CS,'edgecolor','none','verticalalignment','middle',...
    'horizontalalignment','left')     
axis off
subplot(4,3,[8 9 11 12]);
%
g = imagesc(p_vals_mat',[0 1]); axis xy; c=colorbar; hold on; colormap('bone');
colormap( [1 1 1; parula(256)] );caxis( [-0.01 1] );colormap('bone');
xlabel(c,phys_id,'interpreter','latex','fontsize',CS,'rotation',270,'position',[4,0.5,0])

plot(sig_i,sig_j,'r.' ); hold on;
%line([0 length(z_t_given_T)],[0 length(z_t_given_T)],'color','w','linewidth',1.5);
xlabel('Time k (minutes since drug onset) ','fontsize',LS); ylabel('Time j (minutes since drug onset)','fontsize',LS);
%
ax_ticks = time_vec(1):5:time_vec(end);
%
xs = linspace(1,length(time_vec), numel(ax_ticks));
set(gca, 'XTick', xs, 'XTickLabel', ax_ticks);
set(gca,'YTick', xs, 'YTickLabel', ax_ticks);

xl = get(gca,'XLabel');
yl = get(gca, 'YLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', XS)
set(xl, 'FontSize', LS);
yAX = get(gca,'YAxis');
set(yAX, 'FontSize', XS);
set(yl, 'FontSize', LS);

fig_2_id_t = title(fig_2_id,'fontsize',TS,'fontweight','bold');
fig_2_id_t.Position = [xs(1), fig_2_id_t.Position(2:end)];

extensions = {'.pdf','.eps', '.svg'};
for idx = 1:length(extensions)
    ext = extensions{idx};
    %use exportgraphics if you have MATLAB 2020
    %exportgraphics(fig_handle, strcat(figure_params.fn, ext), 'contenttype', 'vector')
    saveas(fig_handle, strcat(figure_params.fn, ext) );
end
%saveas(fig_handle, strcat(figure_params.fn,'.svg'));

end