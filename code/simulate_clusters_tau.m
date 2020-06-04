% Program:  simulate_clusters_tau.m
%
% Summary:  Stochastic simulation of gut bacterial cluster kinetics using
%           the model described in Schlomann and Wiles et al., PNAS (2019).
%           This is a fixed-tau scheme where in each time step the number
%           of each reactions occuring is drawn from a Poisson
%           distribution. There is an addition option to treat growth as a
%           continuous process that is subject to tuneable demographic
%           noise.
%
%           Note: this is code is built off the orginal Growth and
%           Aggregation of Clusters model (github.com/bschloma/gac), which 
%           is now deprecated.
%
% Author:   Brandon Schlomann
%
% Date:     May 2020 --- first written
%
% VCS:      github.com/rplab/cluster_kinetics
%

function [cluster_sizes,total_pop_arr,tvec,num_clusters_arr] = simulate_clusters_tau(growth_rate,aggregation_rate,expulsion_rate,fragmentation_rate,Tmax,n0,max_total_pop,tau,fragmentation_exponent,growth_noise_strength,growth_option)

%% default values for input parameters
% rate of cluster growth
if ~exist('growth_rate','var')||isempty(growth_rate)
    growth_rate = 1;
end

% rate of cluster aggregation
if ~exist('aggregation_rate','var')||isempty(aggregation_rate)
    aggregation_rate = 0.1;
end

% rate of cluster explusion
if ~exist('expulsion_rate','var')||isempty(expulsion_rate)
    expulsion_rate = 0.1;
end

% rate of cluster fragmentation_rate
if ~exist('fragmentation_rate','var')||isempty(fragmentation_rate)
    fragmentation_rate = 0.1;
end

% total simulation time
if ~exist('Tmax','var')||isempty(Tmax)
    Tmax = 24;
end

% initial starting population of single cells
if ~exist('n0','var')||isempty(n0)
    n0 = 10;
end

% carrying capacity
if ~exist('max_total_pop','var')||isempty(max_total_pop)
    max_total_pop = 1e4;
end

% tau
if ~exist('tau','var')||isempty(tau)
    tau = 0.1;
end

% fragmentation exponent: 
% prob_rate_of_frag = fragmentation_rate*(cluster_sizes)^fragmention_exponent*(1-sum(cluster_sizes)/max_total_pop)
if ~exist('fragmentation_exponent','var')||isempty(fragmentation_exponent)
    fragmentation_exponent = 1;
end

% growth_noise_strength = 0 ==> deterministic growth
% growth_noise_strength = 1 ==> exact langvein equation for stochastic
%                               logistic growth, if in 'gaussian' mode
% In 'poisson mode', all values of growth_noise_strength > 0 are equivalent
if ~exist('growth_noise_strength','var')||isempty(growth_noise_strength)
    growth_noise_strength = 0;
end

% discrete ('poisson') or continuous 'gaussian' growth.
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

% main object of the simulation:  array of cluster volumes.
if size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    cluster_sizes = n0;
    num_clusters_arr(1) = numel(n0); 
end


% array for keeping track of total population size
total_pop_arr = zeros(1,num_time_points);
total_pop_arr(1) = sum(cluster_sizes);

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
            if growth_noise_strength > 0
                dBt = sqrt(tau).*randn(1,numel(cluster_sizes));
            else
                dBt = 0;
            end
            
            % Langevin equation for stochastic logistic birth process with
            % tuneable noise strength. 
            cluster_sizes = cluster_sizes + tau.*growth_rate.*cluster_sizes.*(1-sum(cluster_sizes)./max_total_pop) ...
                + growth_noise_strength.*sqrt(growth_rate.*cluster_sizes.*max((1-sum(cluster_sizes)./max_total_pop),0)).*dBt;
            
            cluster_sizes(cluster_sizes<1) = [];
            
        case 'poisson'            
            mean_num_growth_events = max(growth_rate.*tau.*cluster_sizes.*(1-sum(cluster_sizes)./max_total_pop),0);
            
            if growth_noise_strength==0     % determinstic discrete growth
                cluster_sizes = cluster_sizes + round(mean_num_growth_events);
            else
                cluster_sizes = cluster_sizes + poissrnd(mean_num_growth_events,1,numel(cluster_sizes));
            end
    end
    
   % update total population array
   total_pop_arr(s) = sum(cluster_sizes);
   
   % update total number of clumps array
   num_clusters_arr(s) = numel(cluster_sizes);
     
   
   %% fragmentation  
   % do by each cluster
   mean_num_frag_events_for_each_cluster = max(fragmentation_rate.*tau.*cluster_sizes.^fragmentation_exponent.*(1-sum(cluster_sizes)./max_total_pop),0);
   num_frag_events_for_each_cluster  = poissrnd(mean_num_frag_events_for_each_cluster,1,numel(cluster_sizes));
   
   % manually prevent clusters from fragmenting to zero or negative sizes
   num_frag_events_for_each_cluster([cluster_sizes - num_frag_events_for_each_cluster] < 1 ) = 0;
   
   if sum(num_frag_events_for_each_cluster) > 0
       cluster_sizes = cluster_sizes - num_frag_events_for_each_cluster;
       cluster_sizes = [cluster_sizes, ones(1,sum(num_frag_events_for_each_cluster))];       
       num_clusters_arr(s) = numel(cluster_sizes);
    end
    
    %% aggregation
    % determine total number of agg events, then pick clusters randomly
    number_of_possible_agg_reactions = round(.5.*(numel(cluster_sizes).^2 - numel(cluster_sizes)));
    total_prob_rate_of_an_agg_event_happening  = aggregation_rate.*number_of_possible_agg_reactions;  
    num_agg_events_to_happen = poissrnd(total_prob_rate_of_an_agg_event_happening*tau);

    if num_agg_events_to_happen > 0
        % pick pairs of clusters randomly. first generate a 
        %(num clusters x 2) array of random ids. 
%        agg_ids = ceil(numel(cluster_sizes).*rand(num_agg_events_to_happen,2));
        
%         % forbid clusters from aggregating with themselves 
%         agg_ids(diff(agg_ids,[],2)==0,:) = [];
%         
%         % only let clusters aggregate once per time step by identifying
%         % unique pairs
%         [~,unique_rows,~] = unique(agg_ids(:,2));
%         agg_ids = agg_ids(unique_rows,:);
%         agg_ids(ismember(agg_ids(:,2),agg_ids(:,1)),:) = [];
        
        good_ids = 1:numel(cluster_sizes);
        % loop over first clusters and add the size of the second cluster 
        %for a = 1:size(agg_ids,1)
        for a = 1:num_agg_events_to_happen
            
            first_agg_id = good_ids(ceil(numel(good_ids).*rand(1)));
            good_ids_not_first_id = good_ids(good_ids~=first_agg_id);
            
            second_agg_id = good_ids_not_first_id(ceil(numel(good_ids_not_first_id).*rand(1)));
            
%             % if a cluster is chosen to aggregate with itself, draw again
%             ids_are_the_same = agg_ids(a,1)==agg_ids(a,2);
%             %one_cluster_is_nan = sum(isnan(cluster_sizes(agg_ids(a,:))))>0;
%             while ids_are_the_same %|| one_cluster_is_nan
%                 agg_ids(a,:) = good_ids(ceil(numel(good_ids).*rand(1,2)));
%                 ids_are_the_same = agg_ids(a,1)==agg_ids(a,2);
%                 %one_cluster_is_nan = sum(isnan(cluster_sizes(agg_ids(a,:))))>0;
%             end
%             
            % add the second cluster's size to the first cluster
            %cluster_sizes(agg_ids(a,1)) = cluster_sizes(agg_ids(a,1)) + cluster_sizes(agg_ids(a,2)); 
            cluster_sizes(first_agg_id) = cluster_sizes(first_agg_id) + cluster_sizes(second_agg_id);
           
            % mark the second cluster with a NaN for removal later
            % (removing now would require updating agg_ids).
            %cluster_sizes(agg_ids(a,2)) = NaN;
            cluster_sizes(second_agg_id) = NaN;
            
            %good_ids(good_ids==agg_ids(a,2)) = [];
            good_ids(good_ids==second_agg_id) = [];
            
            if numel(good_ids)==1
                break
            end

            % in case the second cluster was chosen to aggregate again with
            % another cluster, re-assign its agg_id to the first cluster's
            %agg_ids(agg_ids==agg_ids(a,2)) = agg_ids(a,1);
            %agg_ids(agg_ids(:,1)==agg_ids(a,2),1) = agg_ids(a,1);
            %agg_ids(agg_ids(:,2)==agg_ids(a,2),2) = agg_ids(a,1);

        end
        
        % remove second cluster from cluster_sizes array
        %cluster_sizes(agg_ids(:,2)) = [];
        cluster_sizes(isnan(cluster_sizes)) = [];
        
        % udpate number of clusters array
        num_clusters_arr(s) = numel(cluster_sizes);
    end
    
    %% expulsion
    % compute total number of expulsion events, then pick clusters randomly
    total_prob_rate_of_explusion_event_happening = expulsion_rate*numel(cluster_sizes);
    num_expulsion_events_to_happen = poissrnd(total_prob_rate_of_explusion_event_happening*tau);

    if num_expulsion_events_to_happen > 0        
        expulsion_ids = unique(ceil(rand(1,num_expulsion_events_to_happen).*numel(cluster_sizes)));       
        cluster_sizes(expulsion_ids) = [];
             
        % update arrays
        num_clusters_arr(s) = numel(cluster_sizes);
        total_pop_arr(s) = sum(cluster_sizes);
    end
    
    if isempty(cluster_sizes)
        return
    end
end
    


   
        
end




 
    
    
    
    
    
    




