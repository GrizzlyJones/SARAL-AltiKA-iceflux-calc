function [ data ] = load40hz( array, filePath, variable, filter )
%LOAD40HZ Summary of this function goes here
%   Detailed explanation goes here

tmpData = ncread(filePath, variable);
            
data = vertcat(array, tmpData(filter(1)));

end

