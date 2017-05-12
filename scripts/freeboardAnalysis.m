function iceSheets = freeboardAnalysis(fluxgate)
%FREEBOARDANALYSIS Splits data between leads in the fluxgate
%   iceSheets = freeboardAnalysis(FLUXGATE) returns a struct with freeboard
%   heights between leads. The freeboard height is calculated as the height
%   of the surface minus the average lead height surrounding the ice. Some
%   sheets of ice are filtered, due to lack of data.
%   The filtration happens if the icesheet has a width of two data points
%   or less, due to insufficient data. It is also filtrated if a NaN value
%   is pressent, due to NaN values being a mask for bodies of water.
%
%   See also INITFLUXGATE

profile = fluxgate.profile;
data = fluxgate.data;

%% Argument checking
if length(data.sla) ~= length(data.class)
    error('Error:\nDimensions don\'' match');
end

%% Generate structs
iceSheets = struct('lon', [], 'lat', [], 'freeboard', [], 'index', []);
iceSheets.index = struct('start', [], 'stop', []);

tmp = struct('heights', [], 'leads', []);

%% init meta data
initial = true;
newFreeboard = false;
index = 1;
indexStop = nan;

%% Lead and ice sheet detection
for i = 1:length(data.sla)
    % if value is a lead and it's a new ice sheet
    if data.class(i) == 4 && ~newFreeboard
        tmpStop = indexStop;
        indexStop = i - 1;
        newFreeboard = true;
        
        % skips the first lead, due to lack of starting lead
        if ~initial
            % filter, see header for more info
            if any(isnan(data.sla(indexStart:indexStop))) || ...
                    ((indexStop - indexStart) < 3)
                continue;
            end
            
            % ice sheet extractin
            iceSheets(index).lon = profile.lon(indexStart:indexStop);
            iceSheets(index).lat = profile.lat(indexStart:indexStop);
            iceSheets(index).index.start = indexStart;
            iceSheets(index).index.stop = indexStop;
            
            tmp(index).heights = data.sla(indexStart:indexStop);
            tmp(index).leads = data.sla(tmpStop:indexStart);
            
            index = index + 1;
        end
    % if value is part of a new ice sheet
    elseif newFreeboard
        indexStart = i;
        newFreeboard = false;
        initial = false;
    end
end
% Includes last leads data
tmp(index).leads = data.sla(tmpStop:indexStart);

%% Freeboard calculations
for i = 1:length(iceSheets)
    avgLead = (mean(tmp(i).leads) + mean(tmp(i+1).leads)) / 2;
    iceSheets(i).freeboard = tmp(i).heights - avgLead;
end
end

