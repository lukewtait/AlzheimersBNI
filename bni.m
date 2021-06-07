function [BNI,P_sz] = bni(net,K,I_0)
% BNI as function of coupling using the theta model
% M.A.Lopes, 2020. Modified by L.Tait 2021. 
% Inputs: 
% - 'net' (connectivity matrix) 
% - 'K'   (global coupling)
% - 'I_0' (distance to SNIC: should be a vector of distances)
%    All I_0 distances will be assumed equal for all nodes.
% 
% See Tait et al. (2021) A Large Scale Brain Network Mechanism for Seizure
% Propensity in AD for details. 

rng('shuffle');

% Fixed parameters:
T=4*10^6;         % # time steps
n_n=5;            % # runs for noise
I_sig=5*1.2*0.1;  % noise level
flag='BNI';           % normalisation flag (BNI=0; NI=1)

N=length(net);    % # nodes
n_I = length(I_0); % # I_0 values

P_sz=zeros(N,n_n,n_I); % initialize P_sz

% Loop over I_0 values
msg = [] ; 
for I_it=1:n_I
    
    % Display output
    fprintf(repmat('\b',1,length(msg))) ; 
    msg = sprintf('Calculating BNI: I_0 = %.4f (%d of %d)',I_0(I_it),I_it,n_I) ; 
    fprintf(msg) ; 
    
    I_0_aux=I_0(I_it)*ones(N,1);
    
    for noise=1:n_n
        [~,P_sz(:,noise,I_it)]=theta_model_S(net,T,K,I_0_aux,I_sig,flag);
    end
    
end

% Calculate BNI, integrate over Psz
mean_P_sz = squeeze(mean(mean(P_sz),2)) ; % average over noise runs and ROIs
BNI = trapz(I_0,mean_P_sz') ; % integrate