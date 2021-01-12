% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts takes the layer connectivity and raw probability of
% impermeability data and calculates the layer prominence metric in 1km
% bins as displayed in figure 2. It is also used to write the final
% combined data file with all the results.
%
% Output:
%   data - Mx7 matrix of doubles
%       data(:,1) - latitude
%       data(:,2) - longitude
%       data(:,3) - layer prominence
%       data(:,4) - probability that layer density exceed 0.81 g/cm^3
%       data(:,5) - layer connectivity score
%       data(:,6) - minimum layer density consistent with radar observations
%       data(:,7) - maximum layer density consistent with radar observations
%
% Data dependencies:
% Probability of impermeability from radar reflectivity inversion - /DemoData/Demo_ImpermProb.mat
% Layer connectivity results - /DemoData/Demo_Connectivity.txt
% Mean climate data - /RefData/Greenland_5km_v1.1.nc
%
% Function dependencies:
% ./ReferenceFunctions/ll2psn
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

% Load the processed connectivity data
tmp =  readmatrix('./DemoData/Demo_Connectivity.txt');
mean_lat = tmp(:,1);
mean_lon = tmp(:,2);
connectivity = tmp(:,3);

% Load the raw probability of impermeability data
load('./DemoData/Demo_ImpermProb.mat');

% Combine results from each radargram on the transect into one continuous
% data structure
data_prob = [];
lat = [];
lon = [];
for k = 1:length(density_prob)
    data_prob = [data_prob density_prob{k}.imperm_prob];
    lat = [lat density_prob{k}.lat]; 
    lon = [lon density_prob{k}.lon]; 
end

% Calculate the approximate depth of the 2012 melt layer in the 2017 data
lat_ref = ncread('Greenland_5km_v1.1.nc', 'lat');
lon_ref = ncread('Greenland_5km_v1.1.nc', 'lon');

line_depth = zeros(size(lat));
for m = 1:length(line_depth)
        opt = abs(lat_ref - lat(m)).^2 + abs(lon_ref - lon(m)).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);
        line_depth(m) = depth(row,col);
end

% Take the max density probability within +/- 8 samples of the expected 2012 layer 
% depth to collapse layer metrics in the vertical 
density_prob = zeros*ones(1, size(data_prob, 2));
for k = 1:length(density_prob)
    start = round((line_depth(k)/0.3) - 8);
    stop = round((line_depth(k)/0.3) + 8);
    if start < 1
        start = 1;
    end
    if stop > 50
        stop = 50;
    end
    
    seg = data_prob(start:stop,k);
    x = find(~isnan(seg));
    if ~isempty(seg(x))
        density_prob(k) = max(seg(x));
    else
        density_prob(k) = NaN;
    end
    if density_prob(k) == 0
        density_prob(k) = NaN;
    end
end

% Find the average probability of impermeability over the same 1 km grid
% bins that were used to calculated the connectivity index
[at, ~, ~, ~] = geodetic_to_along_track(lat, lon, []);
[x,y] = ll2psn(mean_lat, mean_lon);
[x1, y1] = ll2psn(lat, lon);

mean_density_prob = NaN*ones(size(mean_lat));
coverage = zeros(size(mean_lat));
for m = 1:length(mean_lat)
    % Find the observations that fall in the grid bin
    [~, middle] = min(abs(x1 - x(m)).^2 + abs(y1 - y(m)).^2);
    [~, start] = min(abs(at - (at(middle) - 500)));
    [~, stop] = min(abs(at - (at(middle) + 500)));
    tmp1 = density_prob(start:stop);
    q = find(~isnan(tmp1));
    if ~isempty(q) 
        mean_density_prob(m) = mean(tmp1(q));
    end
end

%%%%%%%%%%%%% Calculate the Density Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the raw density data from inverting the radar reflectivities
load('./DemoData/Demo_Density.mat');

% Create a single combined data structure with results from every radargram
% in the flight transect
data_min = [];
data_max = [];
lat = [];
lon = [];
for k = 1:length(density_result)
    data_min = [data_min density_result{k}.inversion_min];
    data_max = [data_max density_result{k}.inversion_max];
    lat = [lat density_result{k}.lat]; 
    lon = [lon density_result{k}.lon]; 
end

% Calculate the approximate depth of the 2012 melt layer of this flight
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

% Find the minimum and maximum consistent density for each trace with +/- 8
% samples of the estimated layer depth
min_density = zeros(1, size(data_min, 2));
max_density = zeros(1, size(data_min, 2));
for k = 1:length(min_density)
    start = round((line_depth(k)/0.3) - 8);
    stop = round((line_depth(k)/0.3) + 8);
    if start < 1
        start = 1;
    end
    if stop > 50
        stop = 50;
    end
    
    seg = data_min(start:stop,k);
    seg2 = data_max(start:stop,k);
    x = find(~isnan(seg));
    if seg(1) ~= 1000 && ~isempty(seg(x))
        min_density(k) = max(seg(x));
        max_density(k) = max(seg2(x));
    elseif seg(1) == 1000
        min_density(k) = 1000;
        max_density(k) = 1000;
    else
        min_density(k) = NaN;
        max_density(k) = NaN;
    end
end

% Average the layer density over 1km bins
[at, ~, ~, ~] = geodetic_to_along_track(lat, lon, []);
[x,y] = ll2psn(mean_lat, mean_lon);
[x1,y1] = ll2psn(lat, lon);

mean_min_density = NaN*ones(size(mean_lat));
mean_max_density = NaN*ones(size(mean_lat));
for m = 1:length(mean_lat)
    % Find the observations that fall in the grid cell
    [~, middle] = min(abs(x1 - x(m)).^2 + abs(y1 - y(m)).^2);
    [~, start] = min(abs(at - (at(middle) - 500)));
    [~, stop] = min(abs(at - (at(middle) + 500)));
    tmp1 = min_density(start:stop);
    tmp2 = max_density(start:stop);
    % Remove NaN or bad roll observations
    q = find(~isnan(tmp1));
    tmp1 = tmp1(q);
    tmp2 = tmp2(q);
    p = find(tmp1 == 1000);
    tmp1(p) = [];
    tmp2(p) = [];
    % Average over grid bin if observations exist
    if ~isempty(tmp1)
        mean_min_density(m) = mean(tmp1);
        mean_max_density(m) = mean(tmp2);
    end
end

%%%%%%%%%%%%%%%%%%%%% Calculate Layer Prominence %%%%%%%%%%%%%%%%%%%%%%%%%%

layer_prom = zeros(size(mean_density_prob));
for k = 1:length(layer_prom)
    layer_prom(k) = (1/2)*(mean_density_prob(k) + connectivity(k));
end


% Remove NaN segments resulting from excluded/lack of data
ind = [];
for k = 1:length(mean_max_density)
    if isnan(mean_max_density(k)) || isnan(mean_density_prob(k))
        ind = [ind k];
    end
end
layer_prom(ind) = [];
mean_density_prob(ind) = [];
mean_min_density(ind) = [];
mean_max_density(ind) = [];
mean_lat(ind) = [];
mean_lon(ind) = [];
connectivity(ind) = [];

%%%%%%%%%%%%%%%%%%% Write Final Data Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = [mean_lat mean_lon layer_prom mean_density_prob connectivity mean_min_density mean_max_density];
writematrix(data, './DerivedData/Demo_FinalResults.txt');