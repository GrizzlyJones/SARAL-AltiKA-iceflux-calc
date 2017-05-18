function [ rslt ] = flattenIcesheet( fluxgate, iceSheet, var )
%FLATTENICESHEET Combines all ice sheets into one array
%   rslt = flattenIceSheet(FLUXGATE, ICEHSEET, VAR) combines alle layers of
%   VAR in ICESHEET and combines it in a single array with length matching
%   FLUXGATE.

rslt = (1:fluxgate.profile.steps) * nan;
for i = 1:length(iceSheet)
    start = iceSheet(i).index.start;
    stop = iceSheet(i).index.stop;
    rslt(start:stop) = iceSheet(i).(var);
end
end

