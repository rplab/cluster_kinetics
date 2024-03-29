% example_gillespie_script.m
%
% Example script showing how to call simulate_clusters_gillespie with a
% range of parameter values.

%% parameters
% time in hours

% logicals for running various parts of the script
l_assemble = true;                              % run the actual simulations
l_save = false;                                 % save the outputs
l_plot = true;                                  % plot the results

maindir = '.';                                  % directory to save outputs
growth_rate = 0.5;
fragmentation_rate = 0.5;   
aggregation_rate = [0, logspace(-5,-3,3)];
expulsion_rate = 0.01;
fragmentation_exponent = 2/3;
expulsion_exponent = 1/3;
aggregation_exponent = 2/3;
K = 1e3;                                        % carrying capacity
timepoints = [24,48,72];                        % time points to save outputs at. 
Tmax = [timepoints(1) diff(timepoints)];        
num_trials = 10;                                % increase this for more replicates
n0 = 10;                                        % initial condition (here, 10 single cells)

% colors for plotting. interpolate between magenta and cyan
reds = linspace(1,0,numel(timepoints));
blues = ones(1,numel(timepoints));
greens = linspace(0,1,numel(timepoints));


%% main loop
if l_assemble
    tic;
    
    % arrays to save outputs
    cluster_sizes_cell = cell(numel(aggregation_rate),num_trials,numel(timepoints));
    reaction_label_cell = cell(numel(aggregation_rate),num_trials,numel(timepoints));
    tvec_cell = cell(numel(aggregation_rate),num_trials,numel(timepoints));
    
    for i = 1:numel(aggregation_rate)
        disp([num2str(i) ' of ' num2str(numel(aggregation_rate))]);
        
        for j = 1:num_trials
            
            for k = 1:numel(timepoints)
                
                if k==1
                    this_n0 = n0;
                else
                    this_n0 = cluster_sizes_cell{i,j,k-1};
                end
                
                % call the simulation function
                [cluster_sizes,~,~,~,reaction_labels] = simulate_clusters_gillespie(growth_rate,aggregation_rate(i),expulsion_rate,fragmentation_rate,Tmax(k),this_n0,aggregation_exponent,fragmentation_exponent,expulsion_exponent,K,false,'','');
                
                cluster_sizes_cell{i,j,k} = cluster_sizes;
                reaction_label_cell{i,j,k} = reaction_labels;
                
                
                
            end
        end
        
        
        
    end
    
    if l_save
        save([maindir filesep 'cluster_sizes_cell'],'cluster_sizes_cell');
        save([maindir filesep 'reaction_label_cell'],'reaction_label_cell');
        save([maindir filesep 'params'],'fragmentation_rate','aggregation_rate','growth_rate','expulsion_rate','fragmentation_exponent','expulsion_exponent','Tmax','timepoints','K')
    end
    runtime = toc;
else 
    % load saved outputs if desired
%     load([maindir filesep 'cluster_sizes_cell'],'cluster_sizes_cell');
%     load([maindir filesep 'params']);       
end

%% plot
if l_plot
    figure('position', [1 480 1440 325]); hold on;
    legendcell = cell(1,numel(timepoints));
    
    % loop over aggregation rate
    for i = 1:numel(aggregation_rate)
        subplot(1,4,i); hold on
        
        % plot a power law guide as a dashed line
        xline = logspace(0,3,5);
        yline = xline.^(-1);
        h = plot(xline,yline,'k--','linewidth',3);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        % loop over time points and plot each curve in a different color
        for t = 1:numel(timepoints)
            
            % color for plotting
            thiscolor = [reds(t),greens(t),blues(t)];

            % pool clusters from all trials
            these_clusters = [];
            for j = 1:size(cluster_sizes_cell,2)
                these_clusters = [these_clusters,cluster_sizes_cell{i,j,t}];
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
        set(gca,'fontsize',24,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0, 1e1, 1e2, 1e3], 'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e-0],'xminortick','off','yminortick','off')
        axis([5e-1 5e3 1e-4 5e0])
        axis square
        title(['\alpha = ' num2str(aggregation_rate(i),2) ' h^{-1}'],'fontsize',24);
        ylabel('{\it{P}}(size > {\it{n}})','fontsize',24)
        xlabel('{\it{n}}, number of cells','fontsize',24)
                
        % add legend to first subplot only
        if i ==1
            legend(legendcell,'location','ne','fontsize',18);
        end
        
  
    end

end

