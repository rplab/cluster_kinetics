% plot_all_cluster_size_cumdist.m
%
% Script for plotting the experimental reverse cumulative distributions
% from the Supplementary .xlsx file.

%% params
% logical for assembling the cumulative distributions from .xlsx file
l_assemble = true;

% strings identifying each strain. matches xlsx file.
strains = {'a01','a02','ent','pls','psd','z36','z20-che','z20-mot'};

% array of numerical ids to include particular strains.
strain_ids = [1,2,3,4,5,6,7,8];

% plot params
marker_cell = {'o','o','o','o','o','o','o'};
color_cell = {[0.8392 0.3373 0.2353],[0.8000 0.3137 0.4078],[0.9294 0.5176 0.2510], ...
    [.72 0 .72], [0.6353 0.8000 0.2431],[0.4 0.4 0.4],  [0 .82 .82], [0.4157 0.3608 0.6196]};
marker_size = 14;
line_width = 4;

%% load data
if ~exist('T','var')
    opts = detectImportOptions('./Supplementary_Data_File.xlsx');
    opts.DataRange = 'A1';
    T = readtable('./Supplementary_Data_File.xlsx',opts,'Sheet','Fig2');
    T.Properties.VariableNames{1} = 'strain';
end

% convert table of sizes to numerical array
cluster_sizes = table2array(T(:,4:end));            % sizes start at column 4 

%% assemble cumulative distributions
if l_assemble
    % cells for pooled and individual distributions
    cum_dist_cell = cell(numel(strain_ids),1);
    indiv_cum_dist_cell = cell(numel(strain_ids),1);
    
    % loop over strains
    for s = 1:numel(strain_ids)
        % append sizes to an empty array
        tmp_sizes = [];
        
        % collect the rows corresponding to this strain
        these_rows = find(strcmp(T.strain,strains{strain_ids(s)})) ;
        
        % allocate a sub-cell to hold the distributions for each sample
        indiv_cum_dist_cell{s} = cell(1,numel(these_rows));
        
        % loop over samples
        for r = 1:numel(these_rows)
            % collect the relevant size data
            these_cluster_sizes = cluster_sizes(these_rows(r),:);
            
            % if there are any NaNs from the table import, remove
            these_cluster_sizes(isnan(these_cluster_sizes)) = [];
            
            % compute cumulative distribution for this sample
            this_indiv_cum_dist = zeros(1,numel(these_cluster_sizes));
            these_cluster_sizes = sort(these_cluster_sizes);
            for n = 1:numel(this_indiv_cum_dist)
                this_indiv_cum_dist(n) = sum(these_cluster_sizes > these_cluster_sizes(n))./numel(these_cluster_sizes);
            end
            
            % store this cumulative distribution along with the sizes used
            % to compute it 
            indiv_cum_dist_cell{s}{r} = [these_cluster_sizes', this_indiv_cum_dist'];
            
            % collect all sizes to be used to compute pooled distribution
            tmp_sizes = [tmp_sizes, these_cluster_sizes];
            
        end
        
        % compute the pooled distribution
        tmp_sizes = sort(tmp_sizes);
        this_cum_dist = zeros(1,numel(tmp_sizes));
        for n = 1:numel(this_cum_dist)
            this_cum_dist(n) = sum(tmp_sizes > tmp_sizes(n))./numel(tmp_sizes);
        end
        
        % store the pooled results
        cum_dist_cell{s} = [tmp_sizes', this_cum_dist'];
        
    end
    
end



%% plot
figure('position',  [ 193    64   747   726]); hold on;
legend_cell = cell(1,numel(strain_ids));
subplot_inds = [1,2,3,4,5,6,7,8];
title_cell = {'{\it{Aeromonas}} ZOR0001', '{\it{Aeromonas}} ZOR0002', '{\it{Enterobacter}} ZOR00014', '{\it{Plesiomonas}} ZOR0011',...
    '{\it{Pseudomonas}} ZWU0006', '{\it{Vibrio}} ZWU0020','{\it{Vibrio}} ZOR0036','{\it{Vibrio}} ZWU0020 {\Delta}che', '{\it{Vibrio}} ZWU0020 {\Delta}mot'};
for s = 1:numel(strain_ids)      
    subplot(3,3,subplot_inds(s)); hold on;
    
    % collect pooled sizes and distribution
    these_sizes = cum_dist_cell{s}(:,1);
    these_cum_dists = cum_dist_cell{s}(:,2);
    
    % plot only unique values to avoid many markers on the plot
    [~,unique_ids,~] = unique(these_sizes);
    these_sizes = these_sizes(unique_ids);
    these_cum_dists = these_cum_dists(unique_ids);
    
%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%% 

    % -1 line
    xline = these_sizes;
    yline = 1e-2.*these_sizes.^(-1);
    plot(xline,yline,'k--','linewidth',3);

    % individual cumulative distributions
    for n = 1:numel(indiv_cum_dist_cell{s})
        indiv_x = indiv_cum_dist_cell{s}{n}(:,1);
        indiv_y = indiv_cum_dist_cell{s}{n}(:,2);
        
        % again, only plot unique values
        [~,unique_ids,~] = unique(indiv_x);
        indiv_x = indiv_x(unique_ids);
        indiv_y = indiv_y(unique_ids);
        
        % plot as transparent markers
        h = scatter(indiv_x,indiv_y,80,color_cell{strain_ids(s)},'filled');
        alpha(h,0.5); 
        
        % connect markers with transparent lines
        p = plot(indiv_x,indiv_y,'linewidth',2,'color',color_cell{strain_ids(s)});
        p.Color(4) = 0.5;
    end
    
    % pooled distributions
    h = scatter(these_sizes,these_cum_dists,160,color_cell{strain_ids(s)},'o','filled','linewidth',1.2,'MarkerEdgeColor','k');
    alpha(h,0.5);
    
    %%%%%%%%%%%%%% style %%%%%%%%%%%%%
    set(gca,'fontsize',16,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0 1e2 1e4],'xminortick','off','yminortick','off')
    axis([0.2 1e4 1e-4 2])
    axis square
    title(title_cell{strain_ids(s)},'fontsize',16);
    
    % axis labels only once
    if s==8
        xlabel('{\it{n }}(number of cells)','fontsize',24)
    end
    
    if s==4
        ylabel('{\it{P}}(size > {\it{n}})','fontsize',24)
    end
     
end


%% plot all pooled distributions on one axis as lines
subplot(3,3,9); hold on;

% loop over strains
for s = 1:numel(strain_ids)
    % collect pooled sizes and distribution
    these_sizes = cum_dist_cell{s}(:,1);
    these_cum_dists = cum_dist_cell{s}(:,2);
    
    % only plot unique values
    [~,unique_ids,~] = unique(these_sizes);
    these_sizes = these_sizes(unique_ids);
    these_cum_dists = these_cum_dists(unique_ids);
    
     % -1 line
    xline = these_sizes;
    yline = 1e-2.*these_sizes.^(-1);
    plot(xline,yline,'k--','linewidth',3);
    
    % cumulative distribution
    plot(these_sizes,these_cum_dists,'-','linewidth',3,'color',color_cell{strain_ids(s)});

end

% style
set(gca,'fontsize',16,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0 1e2 1e4],'xminortick','off','yminortick','off')
axis([0.2 1e4 1e-4 2])
axis square
title('all strains','fontsize',16)