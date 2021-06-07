% Program:  solve_cluster_master_eqn.m
%
% Summary:  Numerically integrate the master equation for the general
%           cluster model on a grid of integer sizes using a simple Euler
%           step. Only outputs final values of cn.
%
% Inputs:   cn0 - (1 x n_max array) initial condition for density of
%               clusters of size n
%           Tmax - time to integrate model up to (hours)
%           dt - numerical timestep
%           beta - fragmentation rate
%           alpha - aggregation rate
%           lambda - expulsion rate
%           nu_F - fragmentation exponent
%           nu_E - expulsion exponent
%           K - carrying capacity
%           nu_A - aggregation exponent
%
% Outpts:   cn - (1 x n_max array) final values of density of cluster of
%               size n
%
% Author:   Brandon Schlomann
%
% VCS:      github.com/rplab/cluster_kinetics


function [cn] = solve_cluster_master_eqn(cn0,Tmax,dt,r,beta,alpha,lambda,nu_F,nu_E,K,nu_A)

% time array
tvec = 0:dt:Tmax;

% initialize density
cn = cn0;

% arrays that are useful for computations
nvec = 1:numel(cn0);                        % n
nvec_m1 = circshift(nvec,1);                % n-1
nvec_p1 = circshift(nvec,-1);               % n+1

for s = 2:numel(tvec)
     
    cn = update_cluster_model_master(cn,dt,r,beta,alpha,lambda,nu_F,nu_E,K,nvec,nvec_m1,nvec_p1,nu_A);
        
    if sum(isnan(cn))>0
        disp(['error: numerical instability'])
        return
    end
    
end



end

function [cn_out] = update_cluster_model_master(cn_prior,dt,r,beta,alpha,lambda,nu_F,nu_E,K,nvec,nvec_m1,nvec_p1,nu_A)

% shifted arrays
cn_m1 = circshift(cn_prior,1);
cn_p1 = circshift(cn_prior,-1);

% total num clusters
M = sum(cn_prior);

% total num cells
N = sum(nvec.*cn_prior);

% convolution term
if alpha > 0
    conv_term = conv((nvec.^nu_A).*cn_prior,(nvec.^nu_A).*cn_prior);
    conv_term = [0, conv_term(1:(numel(cn_prior)-1))];
else
    conv_term = zeros(1,numel(cn_prior));
end

% 1 < n <= nmax
dcn_dt = r.*(1-N./K).*nvec_m1.*cn_m1 + beta.*(1-N./K).*(nvec_p1.^nu_F).*cn_p1 - (r.*(1-N./K).*nvec + beta.*(1-N./K).*nvec.^nu_F + lambda.*nvec.^nu_E).*cn_prior ...
    +0.5.*alpha.*conv_term - alpha.*nvec.^nu_A.*cn_prior.*sum((nvec.^nu_A).*cn_prior);

% n==1
dcn_dt(1) = beta.*(1-N./K).*(nvec_p1(1).^nu_F).*cn_p1(1) - (r.*(1-N./K).*nvec(1) + beta.*(1-N./K).*nvec(1).^nu_F + lambda.*nvec(1).^nu_E).*cn_prior(1) + beta.*(1-N./K).*sum((nvec.^nu_F).*cn_prior) ...
    +0.5.*alpha.*conv_term(1) - alpha.*nvec(1).^nu_A.*cn_prior(1).*sum((nvec.^nu_A).*cn_prior);

% n = nmax
dcn_dt(end) = r.*(1-N./K).*nvec_m1(end).*cn_m1(end)  - (r.*(1-N./K).*nvec(end) + beta.*(1-N./K).*nvec(end).^nu_F + lambda.*nvec(end).^nu_E).*cn_prior(end) ...
    +0.5.*alpha.*conv_term(end) - alpha.*nvec(end).^nu_A.*cn_prior(end).*sum((nvec(end).^nu_A).*cn_prior(end));
    
% integrate
cn_out = cn_prior + dt.*dcn_dt;

end

