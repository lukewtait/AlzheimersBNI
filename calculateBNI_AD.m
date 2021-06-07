function [PLF,BNI,NI,nNI] = calculateBNI_AD(source_eeg)

% Bandpass filter into 4-12 Hz band
source_eeg = bpfilter(source_eeg,4,8) ; 

% Calculate PLF
PLF = plf(source_eeg) ; 

% If BNI or NI aren't requested, return as BNI and NI are computationally
% expensive
if nargout == 1
    return
end

% Default params for BNI/NI calculations
% See Table 2 of Tait et al. (2021) A Large Scale Brain Network Mechanism
% for Seizure Propensity in Alzheimer's disease
I_0 = linspace(-1.7,-0.5,40) ; 
K = 10 ; 

% Calculate BNI
BNI = bni(PLF,K,I_0) ; 

% If NI isn't requested, return as NI calculation is computationally
% expensive
if nargout == 2
    return
end

% Calculate NI
NI = ni(PLF,K,I_0,BNI) ; 

% Calculate nNI
nNI = nni(ni) ; 

end