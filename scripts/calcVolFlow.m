function [ volFlow, flow ] = calcVolFlow( fluxgate, iceSheet, normVelocities )
%CALCFLUX Calculates the volumetric flow rate through the fluxgate
%   volFlow = calcVolFlow(FLUXGATE, ICESHEET, NORMVELOCITIES) calculates 
%   the volumetric flow rate through the FLUXGATE, by integrating the 
%   ICESHEET multiplied with the NORMVELOCITIES. Numerical integration is
%   done by trapezoidal rule.
% 
%   See also INITFLUXGATE, THICKNESS, PROJVELOCITIES

h_i = flattenIcesheet(fluxgate, iceSheet, 'iceThickness');

flow = h_i .* normVelocities;
flow(isnan(flow)) = 0;

volFlow = trapz(fluxgate.profile.cumStep, flow);
end

