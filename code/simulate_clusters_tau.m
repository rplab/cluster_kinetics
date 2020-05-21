% Program:  simulate_clusters_tau.m

function [cluster_sizes,total_pop_arr,tvec,num_clusters_arr] = simulate_clusters_tau(growth_rate,aggregation_rate,expulsion_rate,fragmentation_rate,Tmax,n0,max_total_pop,tau,fragmentation_exponent,sigma,growth_option)

%% default values for input parameters
% rate of cluster growth
if ~exist('growth_rate','var')||isempty(growth_rate)
    growth_rate = 1;
end

% rate of cluster aggregation
if ~exist('aggregation_rate','var')||isempty(aggregation_rate)
    aggregation_rate = 1;
end

% rate of cluster explusion
if ~exist('expulsion_rate','var')||isempty(expulsion_rate)
    expulsion_rate = .1;%.005;
end

% rate of cluster fragmentation_rate
if ~exist('fragmentation_rate','var')||isempty(fragmentation_rate)
    fragmentation_rate = 20;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 72;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end


% carrying capacity
if ~exist('max_total_pop','var')||isempty(max_total_pop)
    max_total_pop = 1e5;
end

% tau
if ~exist('tau','var')||isempty(tau)
    tau = 0.001;
end

% fragmentation_exponent
if ~exist('fragmentation_exponent','var')||isempty(fragmentation_exponent)
    fragmentation_exponent = 0;
end

if ~exist('sigma','var')||isempty(sigma)
    sigma = 0;
end

if ~exist('growth_option','var')||isempty(growth_option)
    growth_option = 'poisson';
end

% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% time
tvec = 0:tau:Tmax;
num_time_points = numel(tvec);

% array for total number of clusters
num_clusters_arr = zeros(1,num_time_points);

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    cluster_sizes = ones(1,round(n0));
    
    % array for keeping track of number of clumps over time
    num_clusters_arr(1) = n0;

elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    cluster_sizes = n0;
    
     % array for keeping track of number of clumps over time
    num_clusters_arr(1) = numel(n0);
    
end


% array for keeping track of total population size
total_pop_arr = zeros(1,num_time_points);
total_pop_arr(1) = sum(n0);

% hour marker
disp_time_marker = 1;
disp_time_increment = 1;


%% fixed parameters

% logical for printing progress to screen
l_print_progress = false;

%% main simulation.  loop over time.
for s = 2:num_time_points 
            
    % print update to console
    if l_print_progress && tvec(s) >= disp_time_marker
        disp(['time = ' num2str(tvec(s),2) ' number of clusters = ' num2str(numel(cluster_sizes))])
        disp_time_marker = disp_time_marker + disp_time_increment;
    end   
           
    %% growth
    switch growth_option
        case 'gaussian'
            if sigma > 0
                % new way
                %[random_numbers, n, lots_of_random_numbers] = select_x_random_numbers(numel(cluster_sizes));
                %dBt = sqrt(tau).*sqrt(2).*erfinv(2.*random_numbers-1);
                
                % old way
                dBt = sqrt(tau).*randn(1,numel(cluster_sizes));
            else
                dBt = 0;
            end
            
            cluster_sizes = cluster_sizes + tau.*growth_rate.*cluster_sizes.*(1-sum(cluster_sizes)./max_total_pop) ...
                + sigma.*sqrt(growth_rate.*cluster_sizes.*max((1-sum(cluster_sizes)./max_total_pop),0)).*dBt;
            
            cluster_sizes(cluster_sizes<1) = [];
            
        case 'poisson'
            
            mean_num_growth_events = max(growth_rate.*tau.*cluster_sizes.^fragmentation_exponent.*(1-sum(cluster_sizes)./max_total_pop),0);
            
            if sigma==0
                cluster_sizes = cluster_sizes + round(mean_num_growth_events);
            else
                cluster_sizes = cluster_sizes + poissrnd(mean_num_growth_events,1,numel(cluster_sizes));
            end
    end
    
   % append total population array
   total_pop_arr(s) = sum(cluster_sizes);
   
   % append total number of clumps array
   num_clusters_arr(s) = numel(cluster_sizes);
     
   
   %% fragmentation
   
   % do by each cluster
   mean_num_frag_events_for_each_cluster = max(fragmentation_rate.*tau.*cluster_sizes.^fragmentation_exponent.*(1-sum(cluster_sizes)./max_total_pop),0);
   num_frag_events_for_each_cluster  = poissrnd(mean_num_frag_events_for_each_cluster,1,numel(cluster_sizes));
   num_frag_events_for_each_cluster([cluster_sizes - num_frag_events_for_each_cluster] < 1 ) = 0;
   
   if sum(num_frag_events_for_each_cluster) > 0
       cluster_sizes = cluster_sizes - num_frag_events_for_each_cluster;
       cluster_sizes = [cluster_sizes, ones(1,sum(num_frag_events_for_each_cluster))];       
       num_clusters_arr(s) = numel(cluster_sizes);
    end
    
    %% aggregation
    % do all in one
    number_of_possible_agg_reactions = round(.5.*(numel(cluster_sizes).^2 - numel(cluster_sizes)));
    total_prob_rate_of_an_agg_event_happening  = aggregation_rate.*number_of_possible_agg_reactions;  
    num_agg_events_to_happen = poissrnd(total_prob_rate_of_an_agg_event_happening*tau);

    if num_agg_events_to_happen > 0
        agg_ids = ceil(numel(cluster_sizes).*rand(num_agg_events_to_happen,2));
        agg_ids(diff(agg_ids,[],2)==0,:) = [];
        [~,unique_rows,~] = unique(agg_ids(:,2));
        agg_ids = agg_ids(unique_rows,:);
        agg_ids(ismember(agg_ids(:,2),agg_ids(:,1)),:) = [];
        
        for a = 1:size(agg_ids,1)
            cluster_sizes(agg_ids(a,1)) = cluster_sizes(agg_ids(a,1)) + cluster_sizes(agg_ids(a,2));
        end
        
        cluster_sizes(agg_ids(:,2)) = [];
        
        num_clusters_arr(s) = numel(cluster_sizes);
    end
    
    %% expulsion
    % 
    total_prob_rate_of_explusion_event_happening = expulsion_rate*numel(cluster_sizes);
    num_expulsion_events_to_happen = poissrnd(total_prob_rate_of_explusion_event_happening*tau);

    if num_expulsion_events_to_happen > 0        
        expulsion_ids = unique(ceil(rand(1,num_expulsion_events_to_happen).*numel(cluster_sizes)));       
        cluster_sizes(expulsion_ids) = [];
             
        num_clusters_arr(s) = numel(cluster_sizes);
        total_pop_arr(s) = sum(cluster_sizes);
    end
    
    if isempty(cluster_sizes)
        return
    end
end
    


   
        
end




 
    
    
    
    
    
    




