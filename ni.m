function NI = ni(net,K,I_0,BNI_initial)
% NI using the theta model
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
flag='NI';           % normalisation flag (BNI=0; NI=1; IP = 1)

net_initial=net; % initial network with all nodes

N=length(net);    % # nodes
n_I = length(I_0); % # I_0 values


P_sz=zeros(N-1,n_n,N,n_I); % initialize


% Loop over nodes
for node_R=1:N
    
    % Remove node from network
    net=net_initial;
    net(node_R,:)=[];
    net(:,node_R)=[];    
    
    % Loop over I_0 values
    for I_it=1:n_I
        
        I_0_aux=I_0(I_it)*ones(N-1,1); 
        
        for noise=1:n_n
            [~,P_sz(:,noise,node_R,I_it)]=theta_model_S(net,T,K,I_0_aux,I_sig,flag);
        end

    end
 
end

% Calculate NI
mean_P_sz = squeeze(mean(mean(P_sz),2)) ; % average over noise runs and ROIs
BNI_final = trapz(I_0,mean_P_sz') ; 

NI = (BNI_initial-BNI_final)./BNI_initial ; 

