function [C_rtrk, pp_start, pp_stop] = waveformAnalysis(waveform, type)
%WAVEFORMANALYSIS Calculate retracked point
%   C_rtrk = waveformAnalysis(WAVEFORM, TYPE) returns retracked bin using
%   either Offset Center of Gravity (OCOG) or Primary Peak Center of
%   Gravity (PP_COG).
%
%   See also OCOG, PRIMARYPEAK.

% Check number of variables and set default retracker to mid power.
if nargin < 2
    type = '';
end

if strcmp(type, 'OCOG')
    % Offset Center of Gravity
    C_rtrk = OCOG(waveform);
elseif strcmp(type, 'PP_COG')
    % Primary peak analysis
    [pp_start, pp_stop] = primaryPeak(waveform);
    C_rtrk = OCOG(waveform(:), pp_start, pp_stop);
else
    % Mid power retracker
    refSign = 0.5;
    [M, I] = max(waveform);
    gateLow = find(waveform(1:I) < M*refSign, 1, 'last');
    C_rtrk = gateLow + (M*refSign - waveform(gateLow)) / (waveform(gateLow+1) - waveform(gateLow));
end
end