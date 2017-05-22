function [ fluxgate ] = interpProfile( fluxgate, gridData, gridVelocity )
%INTERPPROFILE Interpolates from griddata to profile
%   fluxgate = interpProfile(FLUXGATE, GRIDDATA, GRIDVELOCITY) interpolates
%   from a grid to desired fluxgate profile.
%
%   See also INITFLUXGATE

profile = fluxgate.profile;
data = fluxgate.data;

Xq = gridData.Xq;
Yq = griData.Yq;
sla_pp_cog_q = gridData.sla_pp_cog_q;
ssha_q = gridData.ssha_q;
pPq = gridData.pPq;
Wq = gridData.Wq;

%% Interpolates data on path
data.sla = interp2(Xq, Yq, sla_pp_cog_q, profile.lon, profile.lat);
data.pP = interp2(Xq, Yq, pPq, profile.lon, profile.lat);
data.W = interp2(Xq, Yq, Wq, profile.lon, profile.lat);
data.ssha = interp2(Xq, Yq, ssha_q, profile.lon, profile.lat);
data.brng = interp2(Xq, Yq, gridVelocity.brng, profile.lon, profile.lat);
data.magn = interp2(Xq, Yq, gridVelocity.magn, profile.lon, profile.lat);

%% Classifies and masks data
data.sla(~isnan(data.ssha)) = nan;
data.class = zeros(1, profile.steps);
data.class(data.pP > 30 & data.W < 2) = 4;

fluxgate.data = data;
end

