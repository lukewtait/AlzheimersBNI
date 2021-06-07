function source_eeg = bpfilter(source_eeg,flow,fhigh)

% If it is a fieldtrip source structure, convert
if isfield(source_eeg,'pos')
    tmpdat = source_eeg ; 
    source_eeg = struct ; 
    source_eeg.trial{1} = cell2mat(source_eeg.avg.mom(:)) ; 
    source_eeg.time{1} = source.time ; 
    for i = 1:size(source_eeg.trial{1},1)
        source_eeg.label{i} = sprintf('ROI%d',i) ; 
    end
end

% Check fieldtrip is on path
path = fileparts(which('ft_defaults')) ; 
if isempty(path)
    error('Fieldtrip not installed')
end

% Check it is a raw datatype
source_eeg = ft_datatype_raw(source_eeg) ; 

% Bandpass filter
cfg = struct ; 
cfg.bpfilter = 'yes' ; 
cfg.bpfreq = [flow,fhigh] ; 
source_eeg = ft_preprocessing(cfg,source_eeg) ; 