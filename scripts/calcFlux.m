function [ flux ] = calcFlux( fluxgate, iceSheet, normVelocities )
%CALCFLUX Calculates the flux through the fluxgate
%   flux = calcFlux(FLUXGATE, ICESHEET, NORMVELOCITIES) calculates the flux
%   through the FLUXGATE, by integrating the ICESHEET multiplied with the
%   NORMVELOCITIES.
% 
%   See also INITFLUXGATE, THICKNESS, PROJVELOCITIES

h_i = flattenIcesheet(fluxgate, iceSheet, 'iceThickness');

flux = nansum(h_i .* normVelocities .* fluxgate.profile.stepSize);
end

