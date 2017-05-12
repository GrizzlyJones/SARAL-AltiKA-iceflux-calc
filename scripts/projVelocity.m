function [ normVelocity ] = projVelocity( fluxgate )
%PROJVELOCITY Projects the velocities onto the normal of the fluxgate
%   normVelocity = projVelocity(FLUXGATE) project the velocities onto the
%   nomarl of the fluxgate profile. Should be compatible with curved
%   profiles. Based on dot product between vectors.
%   
%   Formula: V = r1 * r2 * cos(theta1 - theta2)
%            V_norm = r * cos(theta - theta_norm)
%
%   See also INITFLUXGATE, INTERPPROFILE, COOR2BRNG

norm = fluxgate.profile.brng - 90;

normVelocity = fluxgate.data.magn .* cos(deg2rad(fluxgate.data.brng - norm));

end

