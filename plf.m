function PLF = plf(source_eeg,freqband)

eeg = source_eeg.trial{1}' ; 
fsample = 1/mean(diff(eeg.time{1})) ; 

%% Correlation network

nCh = size(eeg,2) ;
PLF = zeros(nCh) ;
lags = zeros(nCh) ; 


H = hilbert(eeg);          
phi=angle(H);      
for i = 1:nCh
    dphi = phi-repmat(phi(:,i),1,nCh);
    lags(i,:) = mean(dphi) ; 
    PLF(i,:) = mean(exp(1i*dphi),1);
end
PLF = abs(PLF-eye(nCh));


%% Additional Options

% Zero lag correlations are likely due to volume conduction (Schmidt et al,
% 2014, PLoS Comput Biol 10(11):e1003974). Hence we can set zero lag
% connections to zero. 
thresh = 2*pi*freqband(1)/fsample ; % 1 sample at minimum frequency
PLF(abs(lags)<thresh) = 0 ; 

% Some correlations may be due to indirect connections. To remove these, we
% calculate number of steps in shortest path length between each pair of
% nodes. If the number of steps is > 1, we set the connection to zero. 
numSteps = distance_length(PLF) ; 
PLF = PLF.*(numSteps == 1) ; 


end

%% Nested functions
function numSteps = distance_length(C) ;
    % 
    % Based on function distance_wei.m by:
    %   Mika Rubinov, UNSW/U Cambridge, 2007-2012.
    %   Rick Betzel and Andrea Avena, IU, 2012

    % Distance matrix formed by taking 1/C for each value of C. Hence short
    % correlations have a low distance. 
    xL = C;
    ind = xL~=0;
    xL(ind) = 1./xL(ind);
    xL = xL./max(max(xL));

    % Set up matrices
    n=length(xL);
    D=zeros(n); D(~eye(n))=inf;                 %distance matrix
    numSteps=zeros(n);                                 %number of edges matrix

    % Dijkstra algorithm
    for u=1:n
        S=true(1,n);                            %distance permanence (true is temporary)
        G1=xL;
        V=u;
        while 1
            S(V)=0;                             %distance u->V is now permanent
            G1(:,V)=0;                          %no in-edges as already shortest
            for v=V
                W=find(G1(v,:));                %neighbours of shortest nodes
                [d wi]=min([D(u,W);D(u,v)+G1(v,W)]);
                D(u,W)=d;                       %smallest of old/new path lengths
                ind=W(wi==2);                   %indices of lengthened paths
                numSteps(u,ind)=numSteps(u,v)+1;              %increment no. of edges in lengthened paths
            end

            minD=min(D(u,S));
            if isempty(minD)||isinf(minD),      %isempty: all nodes reached;
                break,                          %isinf: some nodes cannot be reached
            end;

            V=find(D(u,:)==minD);
        end
    end
end