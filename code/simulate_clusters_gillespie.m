% Program:  simulate_clusters_gillespie.m
% 
% Summary:  Cluster kinetics model implemented with
%           Gillespie's direct method.  The main object is an array of
%           cluster sizes.  Clusters are subject to 4 rate processes,
%           aggregation, fragmentation, growth, and dispersal, each with a
%           rate parameter that can depend on cluster size.  Since growth 
%           is the fastest of the rates, and we don't care about
%           demographic stochasticity, it is approximated as a
%           deterministic process.  This means that clusters sizes are no
%           longer integers.
%
%           The simulation starts from a specified initial condition and
%           runs for a specified amount of time.  A subset of the
%           simulation results are output as arrays and an option exists to
%           write every cluster size at every time point to a txt file.
%
%           Based off of gac_gillespie.m. Minor changes, including having
%           carrying capacity regulate fragmentation. Changed order of
%           inputs
%
% Inputs:   growth_rate - (float) growth rate (1/hr)
%           aggregation_rate - (float) aggregation rate (1/hr)
%           expulsion_rate - (float) expulsion rate (1/hr)
%           fragmentation_rate - (float) fragmentation/sprout rate (1/hr)
%           Tmax - (float) simulation time to exit (hr)
%           n0 - (int or 1x[number of clusters] array of floats/ints) specifies initial
%               conditions. if numel(n0)=1, simulation starts with round(n0)
%               monomers.  if numel(n0) > 1, n0 is taken to be the initial
%               cluster size array.
%           nu_A - (float) exponent determining how the aggregation rate
%                   scales with cluster size via (agg rate) = la*(cluster
%                   size)^nu_A.
%           nu_F - (float) exponent determining how the fragmentation rate
%                   scales with cluster size via (frag rate) = la*(cluster
%                   size)^nu_F.
%           max_total_pop - (float, int) carrying capacity
%           lwritetxt - (logical) 0 for not writing to txt file, 1 for yes.
%           txtdir - (str) full path to dir to save txt files
%           txtname - (str) name of txt file
%           nu_E - (float) exponent determining how the expulsion rate
%                   scales with cluster size via (total explusion rate) = explusion_rate*(cluster
%                   size)^nu_E.    
%
% Outputs:  cluster_sizes - (1x(number of clusters) array of floats) sizes of all
%                   clusters at the final time point
%           total_pop_arr - (1x(number of time steps) array of floats)
%                   total population (sum of all cluster sizes) over time
%           tvec - (1x(number of time steps) array of floats) array of
%                   times at which reactions occured
%           num_clumps_arr - (1x(number of time steps)) number of clusters
%                   over time
%
% Author:   Brandon Schlomann
%
% Date:     Summer 2018 - first written
%
% VCS:      github.com/bschloma/gac
%

function [cluster_sizes,total_pop_arr,tvec,num_clumps_arr] = simulate_clusters_gillespie(growth_rate,aggregation_rate,expulsion_rate,fragmentation_rate,Tmax,n0,nu_A,nu_F,nu_E,max_total_pop,lwritetxt,txtdir,txtname)

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

% aggregation exponent
if ~exist('nu_A','var')||isempty(nu_A)
    nu_A = 0;
end

% fragmentation exponent
if ~exist('nu_F','var')||isempty(nu_F)
    nu_F = 0;
end

% carrying capacity
if ~exist('max_total_pop','var')||isempty(max_total_pop)
    max_total_pop = 1e5;
end

% logical for saving to txt file
if ~exist('lwritetxt','var')||isempty(lwritetxt)
    lwritetxt = false;
end

% path to dir for saving txt file
if ~exist('txtdir','var')||isempty(txtdir)
    txtdir = pwd;
end

% name for saving txt file
if ~exist('txtname','var')||isempty(txtname)
    txtname = 'gacout';
end

% exponent for expulsion
if ~exist('nu_E','var')||isempty(nu_E)
    nu_E = 0;
end

% shuffle random number generator
rng('shuffle');


%% intialize some arrays

% main object of the simulation:  array of cluster volumes
if numel(n0)==1
    cluster_sizes = ones(1,round(n0));
    
    % array for keeping track of number of clumps over time
    num_clumps_arr = n0;

elseif size(n0,1) > 1
    disp('gac error:  initial cluster size array n0 must have dim 1xM');
    return
else
    cluster_sizes = n0;
    
     % array for keeping track of number of clumps over time
    num_clumps_arr = numel(n0);
    
end


% array for keeping track of total population size
total_pop_arr = sum(n0);

% pregenerate lots of random numbers for speed.
number_of_random_numbers_to_pre_generate = round(100000);
lots_of_random_numbers = rand(1,number_of_random_numbers_to_pre_generate);

% counter for how many of the pre-generated random numbers have been used.
n = 0;

% time
t = 0;
tvec = t;

% hour marker
disp_time_marker = 1;
disp_time_increment = 1;

% for printing warning about random number usage
l_first_time_running_out_of_randns = 1;

%% fixed parameters
                
% growth timestep
fraction_of_delta_t = 1;
baseline_dt = 0.1;

% include maximum time between reactions. since growth is determinsitc,
% propensity functions may change substantially due to growth (i.e., if
% system starts our with all monomers, probability of fragmentation is
% zero, but as clusters grow, this probability jumps to non-zero. after
% this time, essentially check again to see if any reactions could happen
% sooner than previously predicted.
max_time_between_reactions = min((1/growth_rate)*log(2), Tmax);

% logical for printing progress to screen
l_print_progress = false;

% for writing to txt file, set things up, print a header and time zero data
if lwritetxt
    
    % if desired save directory doesn't exist, make it
    if ~exist(txtdir,'dir')
        mkdir(txtdir);
    end
    
    % handle to txt file, open in append mode
    fid = fopen([txtdir filesep txtname],'a');
    
    % write header
    fprintf(fid,'%s %s\n','time','cluster sizes');
    
    % time zero data to write
    outarr = [t,cluster_sizes];
    
    % format string for data
    format_str = [repmat('%f ',1,numel(outarr)-1), '%f\n'];
    
    % write time zero data
    fprintf(fid,format_str,outarr);
    
    % close the file
    fclose(fid);
    
end
    
%% main simulation.  loop over time.
while t < Tmax 
            
    % print update to console
    if l_print_progress && t >= disp_time_marker
        disp(['time = ' num2str(t,2) ' number of clusters = ' num2str(numel(cluster_sizes))])
        disp_time_marker = disp_time_marker + disp_time_increment;
    end   
          
    % compute probability of expulsion as a function of cluster volume.
    % here, the relationship is linear in linear length dimension.
    lambda_E = expulsion_rate.*(cluster_sizes).^(nu_E);
    
    if sum(lambda_E)==0
        lambda_E = [];
    end
    
    % create an array of ids labelling every possible explusion reaction
    ids_for_expulsion = 1:numel(lambda_E);
    
    % create an accesory array that labels (1) these reactions as "expulsion"
    label_for_lambda_E = 1.*ones(1,numel(lambda_E));
    
    % compute the probability of aggregation in this timestep.
    % A_mat is a lower triangular array that keeps track of the
    % proababilities for every possible aggregation reaction.
    if aggregation_rate > 0
        if nu_A ==0
            lambda_A = aggregation_rate.*(.5.*(numel(cluster_sizes).^2 - numel(cluster_sizes)));
            
            if numel(cluster_sizes) == 1
                lambda_A = [];
            end
            
        else
            A_mat = aggregation_rate.*tril((cluster_sizes'*cluster_sizes).^nu_A,-1);
            
            % extract the non-zero rates into a linear array
            lambda_A = A_mat(A_mat>0);
            
            % keep track of the indices of the non-zero reactions in A_mat's
            % coordinates.  this will be used to identify the clusters in cluster_sizes
            % that are participating in the reaction.
            lin_inds_for_Amat = find(A_mat>0);
        end
    else
        A_mat = 0;
        
        % extract the non-zero rates into a linear array
        lambda_A = A_mat(A_mat>0);
        
        % keep track of the indices of the non-zero reactions in A_mat's
        % coordinates.  this will be used to identify the clusters in cluster_sizes
        % that are participating in the reaction.
        lin_inds_for_Amat = find(A_mat>0);
    end
    
    % double check dimensions
    if size(lambda_A,1) > 1
        lambda_A = lambda_A';
    end
    
    % create labels (2) that keep track of these possible reactions as
    % "aggregation"
    label_for_lambda_A = 2.*ones(1,numel(lambda_A));
    
    % assign ids to each of the possible aggregation reactions.
    ids_for_agg = 1:numel(lambda_A);
        
    % compute probability of fragmentation.  require clusters to be of at least
    % size 2 to fragment. added growth regulation.
    lambda_F = fragmentation_rate.*(cluster_sizes >= 2).*(cluster_sizes).^(nu_F).*(1-sum(cluster_sizes)./max_total_pop);
    
    if sum(lambda_F)==0
        lambda_F = [];
    end
    
    % assign ids to each of the possible fragmentation reactions.
    ids_for_frag = 1:numel(lambda_F);
    
    % assemble all of the reaction ids into a single linear array
    ids_arr = [ids_for_expulsion, ids_for_agg, ids_for_frag];
    
    % create labels (3) that keep track of these possible reactions as
    % "fragmentation" 
    label_for_lambda_F = 3.*ones(1,numel(lambda_F));
    
    % assemble all reaction probability rates into a single linear array
    lambda_arr = [lambda_E, lambda_A, lambda_F];
    
    % compute total proabability rate of a reaction happening
    lambda_total = sum(lambda_arr);
    
    % assemble all labels into a single linear array
    all_labels = [label_for_lambda_E, label_for_lambda_A, label_for_lambda_F]; 
    
    % choose reaction time based on the total probability rate of a
    % reaction happening.
    
    % select a random number from pre generated list.
    % if we've run out, generate some new ones.
    [random_number, n, lots_of_random_numbers] = select_a_random_number();
                        
    delta_t = (1/lambda_total).*log(1/random_number);
    
   
    % choose next reaction using the array of all possible reaction rates
    % BUT if next reaction time is greater than Tmax, don't execute it,
    % only update with growth.  this is denoted with reaction_id = 0.
    if delta_t <= max_time_between_reactions
        
        if t+delta_t <= Tmax
            
            
            % select a random number from pre generated list.
            % if we've run out, generate some new ones.
            [random_number, n, lots_of_random_numbers] = select_a_random_number();
            
            % get the id of the chosen reaction
            reaction_id = find(random_number < cumsum(lambda_arr)./lambda_total,1);
            
            
            
        else
            reaction_id = 0;
            delta_t = Tmax-t;
            
        end
    else
        delta_t = max_time_between_reactions;
        reaction_id = -1;
    end
    
    % if there's a reaction to happen, call the update routine and update
    % the output arrays
    if ~isempty(reaction_id)
             
        [cluster_sizes,tvec,total_pop_arr,num_clumps_arr] = gac_gillespie_update();
        
    end
    
    % update time
    t = t + delta_t;
    
    % if the last reaction was dropped because it was scheduled to happen
    % after Tmax, cull output arrays down to Tmax
    if reaction_id==0
        
        total_pop_arr = total_pop_arr(tvec<=Tmax);
        num_clumps_arr = num_clumps_arr(tvec<=Tmax);
        tvec = tvec(tvec<=Tmax);
        
    end
         
    % if desired, write time and cluster size array to txt file
    if lwritetxt
        
        % open the file in append mode
        fid = fopen([txtdir filesep txtname],'a');
        
        % data to be written
        outarr = [t,cluster_sizes];
        
        % format string for data
        format_str = [repmat('%f ',1,numel(outarr)-1), '%f\n'];
        
        % write the data
        fprintf(fid,format_str,outarr);
        
        % close the file
        fclose(fid);
        
    end
    
    % if the population went extinct (no more clusters), stop the
    % simulation
    if isempty(cluster_sizes)
        return
    end
    
end
    
    % gillespie update function.  nested to inherit random numbers
    function [cluster_sizes_out,tvec_out,tot_pop_out,num_clumps_out] = gac_gillespie_update()
        
        % collect the label that denotes which reaction is happening
        if reaction_id > 0
            this_label = all_labels(reaction_id);
        else
            this_label  = 0;
            
        end
        
        % create and output variable.  necessary to avoid referencing issues
        % with the nested function.
        cluster_sizes_out = cluster_sizes;
        
        %% growth (deterministic) - do this first
        
        % construct a timestep for numerical integration.  can change
        % resolution here for better or worse accuracy.  Base timestep is
        % .01, if time between reactions gets small, use a smaller step.
        dt = min(delta_t.*fraction_of_delta_t,baseline_dt);
        
        % number of steps taken between reactions
        %numsteps = round((delta_t - dt)/dt);
        numsteps = ceil(delta_t/dt);

        
        % update time array.  new variable name is needed to avoid
        % referencing issues with the nested function.
        tvec_out = [tvec, linspace(t+dt, t+delta_t, numsteps)];
            
        % temporary arrays to be updated during numerical integration
        tmp_tot_pop = zeros(1,numsteps);
        tmp_num_clumps = zeros(1,numsteps);
        
        % loop over time and update according to growth equation
        for s = 1:numsteps
        
            % grow all clusters in one deterministic vectorized step
            cluster_sizes_out = cluster_sizes_out + dt.*growth_rate.*cluster_sizes_out.*(1-sum(cluster_sizes_out)./max_total_pop);
            
            % if clusters die out, remove them
            %cluster_sizes_out = cluster_sizes_out(cluster_sizes_out>=1);
            
            % update tmp_totpop
            tmp_tot_pop(s) = sum(cluster_sizes_out);
            
            % update tmp_num_clumps
            tmp_num_clumps(s) = numel(cluster_sizes_out);
        
        end
        
        % append total population array
        tot_pop_out = [total_pop_arr, tmp_tot_pop];
        
        % append total number of clumps array
        num_clumps_out = [num_clumps_arr, tmp_num_clumps];
        
        %% update cluster size array based on which reaction is happening
        switch this_label
            case 1
                % explusion

                % explusion id
                expulsion_id = ids_arr(reaction_id);
                
                % remove clusters that have been expelled from the cluster array.
                cluster_sizes_out(expulsion_id) = [];
                
            case 2
                % aggregation
                              
                if nu_A > 0
                    % find ids of clusters that will aggregate in this timestep
                    agg_id = ids_arr(reaction_id);
                    
                    % get index of this aggregation reaction in A_mat
                    lin_ind_of_this_agg = lin_inds_for_Amat(agg_id);
                    
                
                    % backout the ids of clusters involved in this aggregation
                    % reaction (row + column of A_mat) CAN replace size(A_mat)
                    % with numel(cluster_sizes)^2????
                    [agg_row,agg_col] = ind2sub(size(A_mat),lin_ind_of_this_agg);
                
                elseif nu_A == 0
                    
                    % simply pick two clusters at random
                    
                    l_same_inds = 1;
                    
                    while l_same_inds
                        % select a random number from pre generated list.
                        % if we've run out, generate some new ones.
                        [random_number, n, lots_of_random_numbers] = select_a_random_number();  
                        
                        agg_row = ceil(random_number.*numel(cluster_sizes_out));
                        
                        % select a random number from pre generated list.
                        % if we've run out, generate some new ones.
                        [random_number, n, lots_of_random_numbers] = select_a_random_number();
                                               
                        agg_col = ceil(random_number.*numel(cluster_sizes_out));
                        
                        if agg_col~=agg_row
                            l_same_inds = 0;
                        end
                    end
                
                else
                    disp('error in gac_gillespie: nu_A must be > 0')
                    return
                end
                
                % one of the clusters increases in size due to aggregation
                cluster_sizes_out(agg_row) = cluster_sizes_out(agg_row) + cluster_sizes_out(agg_col);
                
                % the other one is removed from the array
                cluster_sizes_out(agg_col) = [];
                                          
            case 3
                % sprout//fragment
                
                % get id of cluster that will fragment
                frag_id = ids_arr(reaction_id);
                
                % reduce sprouted clusters by 1
                cluster_sizes_out(frag_id) = cluster_sizes_out(frag_id) - 1;
                
                % add this collection of new single cells to cluster array
                cluster_sizes_out = [cluster_sizes_out, 1];
            case 0
                
                % do nothing
                
        end
                
                
        %% final output
        % update population array based on reaction
        tot_pop_out(end) =  sum(cluster_sizes_out);
        
        % update total number of clumps array based on reaction
        num_clumps_out(end) =  numel(cluster_sizes_out); 
        
                  
    end

    function [random_number, n_out, lots_of_random_numbers_out] = select_a_random_number()
        if n + 1 < number_of_random_numbers_to_pre_generate
            random_number = lots_of_random_numbers(n+1);
            n_out = n + 1;
            lots_of_random_numbers_out = lots_of_random_numbers;
        else
            disp('gac: ran out of pre-generated rns, generating new ones')
            
            % generate new ones
            lots_of_random_numbers_out = rand(1,number_of_random_numbers_to_pre_generate);
            
            % collect the one we need
            random_number = lots_of_random_numbers_out(1);
            
            % reset counter
            n_out = 1;
            
        end
        
    end

   
        
end









