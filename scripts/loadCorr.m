function [ data ] = loadCorr( array, filePath, variable, filter )
%LOADCORR Load corrections from file
%   data = loadCorr(ARRAY, FILEPATH, VARIABLE, FILTER) adds corrections
%   (VARIABLE), from FIlEPATH. It then filters the data using FILTER, and
%   adds the data to ARRAY

N = sum(filter(:));
n = sum(filter(1,:));

tmpData = ncread(filePath, variable);

x = linspace(1, n, N);
            
data = vertcat(array, interp1(tmpData(filter(1, :)), x)');

end

