function [ data ] = load40hz( array, filePath, variable, filter )
%LOAD40HZ Summary of this function goes here
%   Detailed explanation goes here

tmpData = ncread(filePath, variable);
            
if ~isempty(strfind(variable, 'waveforms'))
    data = vertcat(array, tmpData(:, filter));
else
data = vertcat(array, tmpData(filter));
end

end

