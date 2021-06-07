% tau_example_script.m
%
% Example of how to call the simulate_clusters_tau.m function with varying parameters.
% Here this function is used to simulate a minimal growth/fragmentation
% process, though it also has functionality for a subset of other cluster
% processes. 



%% params
% time in hours

l_assemble = true;                                      % logical for running simulation
l_save = false;                                         % save outputs
l_plot = true;                                          % plot distribution
maindir = '.';                                          % full path to directory for saving data

r = 0.5;                                                % growth rate
beta = 0.05;                                            % fragmentation rate
timepoints = [6,12,18,24];                              % time points to save output at
Tmax = [timepoints(1) diff(timepoints)];                % max time array for calling main function
n0 = ones(1,10);                                        % initial condition (10 single cells)
tau = 0.01;                                             % algorithm time step
fragmentation_exponent = 1;                             % determines size-dependence of fragmentation (nu_F)
sigma = 1; %                                            % growth noise parameter
growth_option = 'poisson';                              % different types of stochastic growth are available
num_trials = 1;                                         % increase this for more replicates

% unused params
alpha = 0.0;                                            % aggregation rate
lambda = 0.0;                                           % expulsion rate
K = 1/eps;                                              % carrying capacity. approximate exponential growth with large K (1/machine precision)
sig_K = 0.0;                                            % variation in K

% colors for plotting. interpolate between magenta and cyan
reds = linspace(1,0,numel(timepoints));
blues = ones(1,numel(timepoints));
greens = linspace(0,1,numel(timepoints));

%% main loop 

if l_assemble           
    % arrays to save outputs
    cluster_sizes_cell = cell(numel(timepoints),num_trials);
    
    tic;
    % loop over trials
    for m = 1:num_trials
        
        disp([num2str(m) ' of ' num2str(num_trials)]);

        % loop over trials
        for t = 1:numel(Tmax)
            if t > 1
                n0 = cluster_sizes_cell{t-1,m};
            end
                        
            % call the main function
            [cluster_sizes,~,~,~] = simulate_clusters_tau(r,alpha,lambda,beta,Tmax(t),n0,K,tau,fragmentation_exponent,sigma,growth_option);
            
            % collect output into a cell
            cluster_sizes_cell{t,m} = cluster_sizes;
           
        end
        
    end

    if l_save
        save([maindir filesep 'cluster_sizes_cell'],'cluster_sizes_cell');
        save([maindir filesep 'params'],'r','beta','alpha','lambda','timepoints','n0','tau','fragmentation_exponent','sigma','growth_option','K','sig_K');
    end
    
    runtime = toc;
else
      % load previous outputs if desired
%     load([maindir filesep 'cluster_sizes_cell']);
%     load([maindir filesep 'params']);
end

%% plot
if l_plot
    figure('position', [333 280 414 390]); hold on;
    legendcell = cell(1,numel(timepoints));
    
    % plot a power law guide
    xline = logspace(0.5,4.5,5);
    yline = 10.*xline.^(-1);
    h = plot(xline,yline,'k--','linewidth',4);
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    % loop over time points and plot each curve in a different color
    for t = 1:numel(timepoints)
        
        % color for plotting
        thiscolor = [reds(t),greens(t),blues(t)];
        
        % pool clusters from all trials
        these_clusters = [];
        for m = 1:size(cluster_sizes_cell,2)
            these_clusters = [these_clusters,cluster_sizes_cell{t,m}];
        end
        these_clusters = sort(ceil(these_clusters));
        
        % compute reverse cumulative distribution
        this_cum_dist = zeros(1,numel(these_clusters));
        for k = 1:numel(this_cum_dist)
            this_cum_dist(k) = sum(these_clusters > these_clusters(k))./numel(these_clusters);
        end
        
        % plot results
        plot(these_clusters,this_cum_dist,'-','linewidth',3,'color',thiscolor);
        
        % create legend info
        legendcell{t} = ['T = ' num2str(timepoints(t)) ' h'];
        
    end
    
    % style
    set(gca,'fontsize',24,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0, 1e2, 1e4, 1e6], 'ytick',[1e-6 1e-4 1e-2 1e-0],'xminortick','off','yminortick','off')
    axis([5e-1 1e6 1e-6 1e1])
    axis square
    xlabel('{\it{n}}, number of cells','fontsize',24)
    ylabel('{\it{P}}(size > {\it{n}})','fontsize',24)
    title('\nu_F = 1','fontsize',24)
    legend(legendcell,'location','ne','fontsize',16)

end

