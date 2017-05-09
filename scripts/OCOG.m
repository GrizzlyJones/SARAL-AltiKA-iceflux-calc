function [C_rtrk, M, W] = OCOG(waveform, n1, n2)
%OCOG Retrack waveform using Offset Center Of Gravity.
%   C_rtrk = OCOG(WAVEFORM) returns the retracked bin number of
%   the given waveform.
%
%   Source: "Improved sea level determination in the Arctic regions
%            through development of tolerant altimetry retracking"
%           by Maulik Jain

if nargin < 3
    n1 = 1;
    n2 = length(waveform);
end

a = 0;
b = 0;
c = 0;
for i = n1:n2
   a = a + waveform(i)^4;
   b = b + waveform(i)^2;
   c = c + i*waveform(i)^2;
end
M = sqrt(a/b);
W = b^2/a;

COG = c / b;
C_rtrk = COG - W/2;
end