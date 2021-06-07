% plot_all_cluster_size_cumdist_final.m

% logicals
l_assemble = 1;
l_ceil = 1;

% strings identifying each strain. matches xlsx file.
strains = {'a01','a02','ent','pls','psd','z36','z20-che','z20-mot'};

% array of numerical ids to include particular strains.
strain_ids = [1,2,3,4,5,6,7,8];


tot_pop_thresh = 0*100;

% plot params
marker_cell = {'o','o','o','o','o','o','o'};
color_cell = {[0.8392 0.3373 0.2353],[0.8000 0.3137 0.4078],[0.9294 0.5176 0.2510], ...
    [.72 0 .72], [0.6353 0.8000 0.2431],[0.4 0.4 0.4],  [0 .82 .82], [0.4157 0.3608 0.6196]};
marker_size = 14;
line_width = 4;
xbar_width = .1;

% load data
if ~exist('T','var')
    opts = detectImportOptions('/Users/brandonschlomann/Dropbox/acer/Documents/Gutz/clusters/cluster_paper/data/combined_cluster_sizes_markedup_tmp.xlsx');
    opts.DataRange = 'A1';
    %T = readtable('/c/Users/Brandon/Documents/Gutz/clusters/cluster_paper/data/combined_cluster_sizes.xlsx');
    T = readtable('/Users/brandonschlomann/Dropbox/acer/Documents/Gutz/clusters/cluster_paper/data/combined_cluster_sizes_final.xlsx',opts);
    T.Properties.VariableNames{1} = 'strain';
end

cluster_sizes = table2array(T(:,4:end));   

% assemble prob dens
if l_assemble
    cum_dist_cell = cell(numel(strain_ids),1);
    indiv_cum_dist_cell = cell(numel(strain_ids),1);
    
    for s = 1:numel(strain_ids)
     
        tmp_sizes = [];
        
        these_rows = find(strcmp(T.strain,strains{strain_ids(s)}) & T.time == 24) ;
        
        indiv_cum_dist_cell{s} = cell(1,numel(these_rows));
            
        for r = 1:numel(these_rows)
            
            these_cluster_sizes = cluster_sizes(these_rows(r),:);
            these_cluster_sizes(isnan(these_cluster_sizes)) = [];
            
             if l_ceil
                these_cluster_sizes = ceil(these_cluster_sizes);
            end
            
            
           
           this_indiv_cum_dist = zeros(1,numel(these_cluster_sizes));
           these_cluster_sizes = sort(these_cluster_sizes);
           for n = 1:numel(this_indiv_cum_dist)    
               this_indiv_cum_dist(n) = sum(these_cluster_sizes > these_cluster_sizes(n))./numel(these_cluster_sizes);
           end
           
           indiv_cum_dist_cell{s}{r} = [these_cluster_sizes', this_indiv_cum_dist'];
            
            tmp_sizes = [tmp_sizes, these_cluster_sizes];
            
        end
        
        
        tmp_sizes = sort(tmp_sizes);
        this_cum_dist = zeros(1,numel(tmp_sizes));
        for n = 1:numel(this_cum_dist)
            this_cum_dist(n) = sum(tmp_sizes > tmp_sizes(n))./numel(tmp_sizes);
        end
        
        cum_dist_cell{s} = [tmp_sizes', this_cum_dist'];
        
    end
    
end



%% plot

figure('position',  [ 193    64   747   726]); hold on;
legend_cell = cell(1,numel(strain_ids));
subplot_inds = [1,2,3,4,5,6,7,8];
title_cell = {'{\it{Aeromonas}} ZOR0001', '{\it{Aeromonas}} ZOR0002', '{\it{Enterobacter}} ZOR00014', '{\it{Plesiomonas}} ZOR0011',...
    '{\it{Pseudomonas}} ZWU0006', '{\it{Vibrio}} ZWU0020','{\it{Vibrio}} ZOR0036','{\it{Vibrio}} ZWU0020 {\Delta}che', '{\it{Vibrio}} ZWU0020 {\Delta}mot'};
integrated_diffs = [];
for s = 1:numel(strain_ids)
        
    subplot(3,3,subplot_inds(s)); hold on;
    
    these_sizes = cum_dist_cell{s}(:,1);
    these_cum_dists = cum_dist_cell{s}(:,2);
    
    [~,unique_ids,~] = unique(these_sizes);
    these_sizes = these_sizes(unique_ids);
    these_cum_dists = these_cum_dists(unique_ids);
    
%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%% 

    % -1 line
    xline = these_sizes;
    yline = 1e-2.*these_sizes.^(-1);
    plot(xline,yline,'k--','linewidth',3);

    % indiv cum dist
    
    for n = 1:numel(indiv_cum_dist_cell{s})
        
        indiv_x = indiv_cum_dist_cell{s}{n}(:,1);
        indiv_y = indiv_cum_dist_cell{s}{n}(:,2);
        
        [~,unique_ids,~] = unique(indiv_x);
        indiv_x = indiv_x(unique_ids);
        indiv_y = indiv_y(unique_ids);
        
        h = scatter(indiv_x,indiv_y,80,color_cell{strain_ids(s)},'filled');
        alpha(h,0.5); 
        
        p = plot(indiv_x,indiv_y,'linewidth',2,'color',color_cell{strain_ids(s)});
        p.Color(4) = 0.5;
    end
    
    % cum dist
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


%% one axis

subplot(3,3,9); hold on;

for s = 1:numel(strain_ids)
    
    these_sizes = cum_dist_cell{s}(:,1);
    these_cum_dists = cum_dist_cell{s}(:,2);
    
    [~,unique_ids,~] = unique(these_sizes);
    these_sizes = these_sizes(unique_ids);
    these_cum_dists = these_cum_dists(unique_ids);
    
     % -1 line
    xline = these_sizes;
    yline = 1e-2.*these_sizes.^(-1);
    plot(xline,yline,'k--','linewidth',3);
    
    % cum dist
    plot(these_sizes,these_cum_dists,'-','linewidth',3,'color',color_cell{strain_ids(s)});

    
end

set(gca,'fontsize',16,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0 1e2 1e4],'xminortick','off','yminortick','off')
axis([0.2 1e4 1e-4 2])
axis square

title('all strains','fontsize',16)