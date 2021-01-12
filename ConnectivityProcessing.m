% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts take the raw layer detections from the radar data,
% calculates the connectivity metric, and averages it over 1 km bins to
% create the final connectivity dataset shown in Figure 2.
%
% Outputs:
%   data - Mx3 matrix of doubles with location and connectivity scores for
%   each 1 km grid bin
%       data(:,1) - latitude
%       data(:,2) - longitude
%       data(:,3) - connectivity score
%
% Data dependencies:
% Raw layer detections from radar - /DemoData/Demo_Layer.mat
% Mean climate data - /RefData/Greenland_5km_v1.1.nc
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Calculate the approximate expected depth of the 2012 summer surface using
% the steady state Herron & Langway model - to be used later for rough
% localization of the melt layer in the radar data

rho_ice = 0.917;          % density of solid ice in g/cm^3
R = 8.314;                % gas constant
rho1 = 0.33:0.001:0.55;   % density axis

% Load mean annual accumulation rates and surface temperatures
accum = ncread('./RefData/Greenland_5km_v1.1.nc', 'presprcp');
temp = ncread('./RefData/Greenland_5km_v1.1.nc', 'surftemp');

% Convert accumulation from mWe to mIe
accum = accum*(997/917);

depth = zeros(size(accum));
for m = 1:size(accum,1)
    for n = 1:size(accum,2)
        
        A = accum(m,n);
        T = 273.15 + temp(m,n);
        
        % Herron & Langway steady state solution for depth and age over the
        % firn in the 1st compaction regime
        k0 = 11*exp(-10160/(R*T));
        h1 = (1./(rho_ice.*k0)).*(log(rho1./(rho_ice-rho1)) - log(rho1(1)./(rho_ice - rho1(1))));
        t1 = (1/(k0*A)).*log((rho_ice - rho1(1))./(rho_ice - rho1));
        
        % Calculate the depth of the 2012 horizon in 2017 (~4.66 years
        % after initial deposition)
        [~, ind] = min(abs(t1 - 4.66));
        depth(m,n) = h1(ind);
    end
end

% Load the raw radar ice layer detections and sum up the number of
% horizontally adjacent detections for each pixel 
load('./DemoData/Demo_Layer_Connectivity.mat');

% Clip the data file to the region around the surface (may need to manually
% tweak this depending on the flight altitude variations)
upper = 420;
lower = 1700;

% Combine results from each radargram on the transect into one continuous
% data structure
data = [];
lat = [];
lon = [];
surf_ind = [];
for k = 1:length(results)
    data = [data results{k}.data2(upper:lower,:)];
    lat = [lat results{k}.lat];
    lon= [lon results{k}.lon];
    surf_ind = [surf_ind surf{k}];
end

hc = zeros(size(data));      % measured connectivity
pc = NaN*ones(size(data));   % connectivity of a perfectly continuous layer
for m = 2:size(hc,1)-1
    for n = 2:size(hc,2)-1
        
        % If data in a pixel has been rejected, than the connectivity is NaN
        if isnan(data(m,n))
            hc(m,n) = NaN;
        end
        
        % Set ideal score to 3 if a layer is detected
        if ~isnan(data(m,n))
            pc(m,n) = 3;
        end
        % If the the left or right set of pixels is NaN (due to data
        % rejected for aircraft roll), reduce the ideal score accordingly
        if isnan(data(m,n-1)) && data(m,n) ~= -200
            pc(m,n) = pc(m,n) - 1;
        end
        if isnan(data(m,n+1)) && data(m,n) ~= -200
            pc(m,n) = pc(m,n) - 1;
        end
        
        % If a layer has been detected, then sum the number of horizontally
        % adjacent detections
        if data(m,n) ~= -200 && ~isnan(data(m,n)) && data(m,n) ~= 0
            hc(m,n) = 1;
            if data(m,n+1) ~= -200 && ~isnan(data(m,n+1)) && data(m,n+1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
            if data(m,n-1) ~= -200 && ~isnan(data(m,n-1)) && data(m,n-1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
            if data(m-1,n+1) ~= -200 && ~isnan(data(m-1,n+1)) && data(m-1,n+1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
            if data(m+1,n+1) ~= -200 && ~isnan(data(m+1,n+1)) && data(m+1,n+1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
            if data(m-1,n-1) ~= -200 && ~isnan(data(m-1,n-1)) && data(m-1,n-1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
            if data(m+1,n-1) ~= -200 && ~isnan(data(m+1,n-1)) && data(m+1,n-1) ~= 0
                hc(m,n) = hc(m,n) + 1;
            end
        end
    end
end

hc_total_2013 = hc;
perfect_score_2013 = pc;

% Calculate the approximate depth of the 2012 horizon over the flight
% transect
lat_ref = ncread('Greenland_5km_v1.1.nc', 'lat');
lon_ref = ncread('Greenland_5km_v1.1.nc', 'lon');

line_depth = zeros(size(lat));
for m = 1:length(line_depth)
        opt = abs(lat_ref - lat(m)).^2 + abs(lon_ref - lon(m)).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);
        line_depth(m) = depth(row,col);
end

% Calculate the approximate data index of the layer in each trace
surf_ind = round(surf_ind - upper + 1 + line_depth/0.4 + 2);

% Calculate the 1km bin average of the raw connectivity scores

% Bin the flight transect in to 1km long grid cells
[at, ~, ~, ~] = geodetic_to_along_track(lat, lon, []);
x = 0:1000:at(end);

depth_hc = NaN*ones(1, length(x));
mean_lat = zeros(size(x));
mean_lon = zeros(size(x));
for m = 1:length(x) - 1
    % Find the observations that fall within that grid cell
    [~, start] = min(abs(at - x(m)));
    [~, stop] = min(abs(at - x(m + 1)));
    % Find the coordinates of the grid cell center
    mean_lat(m) = mean(lat(start:stop));
    mean_lon(m) = mean(lon(start:stop));
    % As long as there are at least 10 observations in the grid cell,
    % calculate the average connectivity score at every depth
    if stop > start + 10 
        ref = zeros(1, stop-start+1);
        meas = zeros(1, stop-start+1);
        for k = start:stop
            % Total score if there was a perfectly connected layer at that
            % depth
            ref(k-start+1) = max(perfect_score_2013(surf_ind(k)-4:surf_ind(k)+4,k));
            % Measured value normalized by the perfect score
            meas(k-start+1) = max(hc_total_2013(surf_ind(k)-4:surf_ind(k)+4,k));
        end
        % Average normalized connectivity score in the grid bin
        depth_hc(m) = nanmean(meas./ref);
    end
end

connectivity = depth_hc;

% Remove any segments where the connectivity is NaN (no observations due to
% aircraft roll or lack of layer detections
ind = [];
for k = 1:length(connectivity)
    if isnan(connectivity(k))       
        ind = [ind k];
    end
end
connectivity(ind) = [];
mean_lat(ind) = [];
mean_lon(ind) = [];

figure;
plot(connectivity);

% Write final data to text file
data = [mean_lat' mean_lon' connectivity'];
writematrix(data, './DerivedData/Demo_Connectivity.txt');

