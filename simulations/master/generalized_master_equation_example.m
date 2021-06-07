% generalized_master_equation_example.m
%
% example of how to call the generalized master equation solver with variable inputs.

%% params
% time in hours

% logicals for running various parts of the script
l_assemble = true;                      % run the numerical integration
l_save = false;                         % save the output
l_plot = true;                          % plot the output

maindir = '.';                          % directory for saving output data here
r = 0.5;                                % growth rate
beta = 0.5;                             % fragmentation rate
alpha = [10^(-2), 10^(-2.5), 10^(-3)];  % aggregation rate, one per nu_A
nu_A = [0,1/3,2/3];                     % aggregation exponent
lambda = 0.01;                          % expulsion rate
nu_F = 2/3;                             % fragmentation exponent
nu_E = 1/3;                             % expulsion exponent
dt = 0.0001;                            % time step
K = 1e2;                                % carrying capacity
Tmax = 24;                              % simulation time

nmax = ceil(K);                         % max size for master equation integration


% colors for plotting
reds = linspace(1,0,numel(nu_F)*numel(alpha));
blues = ones(1,numel(nu_F)*numel(alpha));
greens = linspace(0,1,numel(nu_F)*numel(alpha));


%% solve master equation
if l_assemble
    
    % initial condition: 10 single cells
    cn0 = zeros(1,nmax);
    cn0(1) = 10;
    
    % store output in a cell
    cn_cell = cell(numel(nu_A),1);
    
    for i = 1:numel(nu_A)
                   
        disp([num2str(i) ' of ' num2str(numel(nu_A))]);
        
        [cn_cell{i}] = solve_cluster_model_master_eqn_final_only_nuA(cn0,Tmax,dt,r,beta,alpha(i),lambda,nu_F,nu_E,K,nu_A(i));
                  
    end
    
    if l_save
        save([maindir filesep 'cn_cell'],'cn_cell');
        save([maindir filesep 'params'],'beta','alpha','r','lambda','nu_F','nu_E','nu_A','dt','Tmax','timepoints','K','nmax')
    end
else
      % load saved outputs if desired
%     load([maindir filesep 'cn_cell.mat'])
%     load([maindir filesep 'params.mat'])
end

%% plot
if l_plot
    figure('position', [117 483 964 287]); hold on;
    nvec = 1:nmax;
    for i = 1:numel(nu_A)
        subplot(1,3,i); hold on
        thiscolor = [reds(i),greens(i),blues(i)];
        
        this_cn = cn_cell{i};
        
        % normalize to get probability, than sum to get cumulative dist
        this_prob_dens = this_cn./sum(this_cn);
        this_cum_dist = 1-cumsum(this_prob_dens);
        
        plot(nvec,this_cum_dist,'linewidth',3,'color',thiscolor);
        set(gca,'fontsize',24,'linewidth',4,'xscale','log','yscale','log','xtick',[1e0, 1e1, 1e2, 1e3], 'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e-0],'xminortick','off','yminortick','off')
        axis([5e-1 5e3 5e-5 5e0])
        axis square
        title(['\nu_A=' num2str(nu_A(i),2)],'fontsize',24);
        ylabel('{\it{P}}(size > {\it{n}})','fontsize',24);
        xlabel('{\it{n}}, number of cells','fontsize',24)
        
    end
end


