function [pp_start, pp_stop, Th_start, Th_stop, d_1] = primaryPeak(waveform)
%primaryPeak Calculate primary peak 
%   [pp_start, pp_stop, Th_start, Th_stop, d_1] = PRIMARYPEAK(WAVEFORM)
%   Returns the start and stop bin of the primary peak (pp_start, pp_stop)
%   Optional output include the start and stop threshold (Th_start, Th_stop
%   respectively) and d_1.
%
%   Source: "Improved sea level determination in the Arctic regions
%            through development of tolerant altimetry retracking"
%           by Maulik Jain

pp_start = -1;
pp_stop = -1;

N = length(waveform);

d_2 = zeros(N-2,1);
d_1 = zeros(N-2,1);

for i = 1:(N-2);
    d_2(i) = waveform(i + 2) - waveform(i);
    d_1(i) = waveform(i + 1) - waveform(i);
end

Th_start = sqrt(((N-2)*sum(d_2.^2) - sum(d_2)^2)/((N-2)*(N-3)));
Th_stop = sqrt(((N-1)*sum(d_1.^2) - sum(d_1)^2)/((N-1)*(N-2)));

for i = 1:(N-2);
    if d_1(i) > Th_start && pp_start == -1
        pp_start = i - 2;
    end
    if d_1(i) < Th_stop && pp_start ~= -1
        pp_stop = i + 2;
        break;
    end
end
end