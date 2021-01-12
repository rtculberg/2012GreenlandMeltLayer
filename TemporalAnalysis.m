% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts calculates the layer connectivity and density along the
% overlapping 2013 transects and plots the data in the format used for
% Figure 4. It relies on the same base code as the ConnectivityProcessing.m
% and PostProcessing.m scripts.
%
% Output:
%   Transects B and C in Figure 4 of main text
%
% Data dependencies:
% Raw layer connectivity results (2013) - /DemoData/Demo_Connectivity_2013.txt
% Raw layer density results (2013) - /DemoData/Demo_Density_2013.txt
% Final results on 2017 transect - /DemoData/Demo_FinalResults.txt
% Mean climate data - /RefData/Greenland_5km_v1.1.nc
% MAR-derived melt production - /RefData/TotalMelt_2013-2017.mat
%
% Function dependencies:
% ./ReferenceFunctions/ll2psn
% ------------------------------------------------------------------------
tic
clear;

addpath(genpath('./ReferenceFunctions/'));

% Calculate the approximate expected depth of the 2012 summer surface using
% the steady state Herron & Langway model - to be used later for rough
% localization of the melt layer in the radar data

rho_ice = 0.917;          % density of solid ice in g/cm^3
R = 8.314;                % gas constant
rho1 = 0.33:0.001:0.55;   % density axis

% Load mean annual accumulation rates and surface temperatures
accum = ncread('Greenland_5km_v1.1.nc', 'presprcp');
temp = ncread('Greenland_5km_v1.1.nc', 'surftemp');

% Convert accumulation from mWe to mIe
accum = accum*(997/917);

depth_2013 = zeros(size(accum));
for m = 1:size(accum,1)
    for n = 1:size(accum,2)
        
        A = accum(m,n);
        T = 273.15 + temp(m,n);
        
        % Herron & Langway steady state solution for depth and age over the
        % firn in the 1st compaction regime        
        k0 = 11*exp(-10160/(R*T));
        h1 = (1./(rho_ice.*k0)).*(log(rho1./(rho_ice-rho1)) - log(rho1(1)./(rho_ice - rho1(1))));
        t1 = (1/(k0*A)).*log((rho_ice - rho1(1))./(rho_ice - rho1));
        
        % Calculate the depth of the 2012 horizon in spring 2013 (~0.66 years
        % after initial deposition)
        [~, ind] = min(abs(t1 - 0.66)); 
        depth_2013(m,n) = h1(ind);
    end
end


% Load the final results from the 2017 radar data on this transect
data = readmatrix('./DemoData/Demo_FinalResults.txt');

lat_2017 = data(:,1);
lon_2017 = data(:,2);
mean_max_density_2017 = data(:,7);
mean_min_density_2017 = data(:,6);
connectivity_2017 = data(:,5);

% Load the raw density inversion data for the 2013 flight on this transect
load('./DemoData/Demo_Density_2013.mat');

% Combine results from each radargram on the transect into one continuous
% data structure
data_min_2013 = [];
data_max_2013 = [];
bg_2013 = [];
lat_2013 = [];
lon_2013 = [];
for k = 1:length(density_result)
    data_min_2013 = [data_min_2013 density_result{k}.inversion_min];
    data_max_2013 = [data_max_2013 density_result{k}.inversion_max];
    bg_2013 = [bg_2013 density_result{k}.base_density];
    lat_2013 = [lat_2013 density_result{k}.lat]; 
    lon_2013 = [lon_2013 density_result{k}.lon]; 
end

% Calculate the approximate depth of the 2012 melt layer in the 2013 data
lat_ref = ncread('Greenland_5km_v1.1.nc', 'lat');
lon_ref = ncread('Greenland_5km_v1.1.nc', 'lon');

line_depth_2013 = zeros(size(lat_2013));
for m = 1:length(line_depth_2013)
        opt = abs(lat_ref - lat_2013(m)).^2 + abs(lon_ref - lon_2013(m)).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);
        line_depth_2013(m) = depth_2013(row,col);
end

% Take the max density within +/- 8 samples of the expected 2012 layer 
% depth to collapse layer metrics in the vertical 
min_density_2013 = zeros(1, size(data_min_2013, 2));
max_density_2013 = zeros(1, size(data_min_2013, 2));
for k = 1:length(min_density_2013)
    start = round((line_depth_2013(k)/0.3) - 8);
    stop = round((line_depth_2013(k)/0.3) + 8);
    if start < 1
        start = 1;
    end
    if stop > 50
        stop = 50;
    end
    
    seg = data_min_2013(start:stop,k);
    seg2 = data_max_2013(start:stop,k);
    x = find(~isnan(seg));
    if seg(1) ~= 1000 && ~isempty(seg(x))
        min_density_2013(k) = max(seg(x));
        max_density_2013(k) = max(seg2(x));
    elseif seg(1) == 1000
        min_density_2013(k) = 1000;
        max_density_2013(k) = 1000;
    else
        min_density_2013(k) = NaN;
        max_density_2013(k) = NaN;
    end
end

% Build a common along-track axis based on the 2017 flight transect
[x_2017, y_2017] = ll2psn(lat_2017, lon_2017);

% Clip transect to the overlapping segments of interest
lat11 = 69.117861;    % start point of overlap #1
lon11 = -33.611442;
lat12 = 67.809925;    % end point of overlap #1
lon12 = -44.239811;

lat21 = 66.359673;    % start point of overlap #2
lon21 = -45.665488;
lat22 = 66.313979;    % start point of overlap #2
lon22 = -41.331289;

% Convert to polar stereographic coordinates
[x11, y11] = ll2psn(lat11, lon11);
[x12, y12] = ll2psn(lat12, lon12);
[x21, y21] = ll2psn(lat21, lon21);
[x22, y22] = ll2psn(lat22, lon22);

% Find index of start and stop points
[~, start1] = min(abs(x_2017 - x11).^2 + abs(y_2017 - y11).^2);
[~, stop1] = min(abs(x_2017 - x12).^2 + abs(y_2017 - y12).^2);
[~, start2] = min(abs(x_2017 - x21).^2 + abs(y_2017 - y21).^2);
[~, stop2] = min(abs(x_2017 - x22).^2 + abs(y_2017 - y22).^2);

% Generate new along-track axis for each overlapping segments
at_axis1 = [x_2017(start1:stop1) y_2017(start1:stop1)];
at_axis2 = [x_2017(start2:stop2) y_2017(start2:stop2)];

% Clip and combine 2017 data along these two segments of interest
mean_min_density_2017 = vertcat(mean_min_density_2017(start1:stop1),mean_min_density_2017(start2:stop2));
mean_max_density_2017 = vertcat(mean_max_density_2017(start1:stop1),mean_max_density_2017(start2:stop2));
connectivity_2017 = vertcat(connectivity_2017(start1:stop1),connectivity_2017(start2:stop2));

% Final along-track axis
at_axis = vertcat(at_axis1, at_axis2);

% Bin the 2013 density data to the same 1km bins as the 2017 data
[x_2013, y_2013] = ll2psn(lat_2013, lon_2013);

mean_lat = zeros(length(at_axis),1);
mean_lon = zeros(length(at_axis),1);
mean_min_density_2013 = NaN*ones(length(at_axis),1);
mean_max_density_2013 = NaN*ones(length(at_axis),1);
for m = 1:length(at_axis) - 1
    % Find observations in the grid bin
    [~, start13] = min(abs(x_2013 - at_axis(m,1)).^2 + abs(y_2013 - at_axis(m,2)).^2);
    [~, stop13] = min(abs(x_2013 - at_axis(m+1,1)).^2 + abs(y_2013 - at_axis(m+1,2)).^2);
    mean_lat(m) = mean(lat_2013(start13:stop13));
    mean_lon(m) = mean(lon_2013(start13:stop13));
    tmp3 = min_density_2013(start13:stop13);
    tmp4 = max_density_2013(start13:stop13);
    % Excluded NaN or excess roll observations
    q2 = find(~isnan(tmp3));
    tmp3 = tmp3(q2);
    tmp4 = tmp4(q2);
    p2 = find(tmp3 == 1000);
    tmp3(p2) = [];
    tmp4(p2) = [];
    % Average over grid bin
    if ~isempty(tmp3)
        mean_min_density_2013(m) = mean(tmp3);
        mean_max_density_2013(m) = mean(tmp4);
    end
end

% New along-track axis
[at_avg, ~, ~, ~] = geodetic_to_along_track(mean_lat, mean_lon, []);

% Load the raw layer detections for the 2013 transect
load('./DemoData/Demo_Layer_Connectivity_2013.mat');

% Clip data to the area around the surface and layer (this may need to be
% manually tweaked depending on the flight altitude profile)
upper = 1250;
lower = 2000;

% Combine results from each radargram on the transect into one continuous
% data structure
data_2013 = [];
lat_2013 = [];
lon_2013 = [];
surf_ind = [];
for k = 1:length(results)
    data_2013 = [data_2013 results{k}.data2(upper:lower, :)];
    lat_2013 = [lat_2013 results{k}.lat];
    lon_2013= [lon_2013 results{k}.lon];
    surf_ind = [surf_ind surf{k}];
end

hc_2013 = zeros(size(data_2013));  % measured connectivity
pc_2013 = zeros(size(data_2013));  % connectivity of a perfectly continuous layer
for m = 2:size(hc_2013,1)-1
    for n = 2:size(hc_2013,2)-1
        
        % If data in a pixel has been rejected, than the connectivity is NaN
        if isnan(data_2013(m,n))
            hc_2013(m,n) = NaN;
        end
        
        % Set ideal score to 3 if a layer is detected
        if ~isnan(data_2013(m,n))
            pc_2013(m,n) = pc_2013(m,n) + 3;
        end
        % If the the left or right set of pixels is NaN (due to data
        % rejected for aircraft roll), reduce the ideal score accordingly
        if isnan(data_2013(m,n-1))
            pc_2013(m,n) = pc_2013(m,n) - 1;
        end
        if isnan(data_2013(m,n+1))
            pc_2013(m,n) = pc_2013(m,n) - 1;
        end
        
        % If a layer has been detected, then sum the number of horizontally
        % adjacent detections        
        if data_2013(m,n) ~= -200 && ~isnan(data_2013(m,n)) && data_2013(m,n) ~= 0
            hc_2013(m,n) = hc_2013(m,n) + 1;
            if data_2013(m,n+1) ~= -200 && ~isnan(data_2013(m,n+1)) && data_2013(m,n+1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
            if data_2013(m,n-1) ~= -200 && ~isnan(data_2013(m,n-1)) && data_2013(m,n-1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
            if data_2013(m-1,n+1) ~= -200 && ~isnan(data_2013(m-1,n+1)) && data_2013(m-1,n+1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
            if data_2013(m+1,n+1) ~= -200 && ~isnan(data_2013(m+1,n+1)) && data_2013(m+1,n+1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
            if data_2013(m-1,n-1) ~= -200 && ~isnan(data_2013(m-1,n-1)) && data_2013(m-1,n-1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
            if data_2013(m+1,n-1) ~= -200 && ~isnan(data_2013(m+1,n-1)) && data_2013(m+1,n-1) ~= 0
                hc_2013(m,n) = hc_2013(m,n) + 1;
            end
        end
    end
end

hc_total_2013 = hc_2013;
perfect_score_2013 = pc_2013;

% Calculate the approximate depth of the 2012 layer in 2013 on this
% transect
lat = ncread('Greenland_5km_v1.1.nc', 'lat');
lon = ncread('Greenland_5km_v1.1.nc', 'lon');

line_depth_2013 = zeros(size(lat_2013));
for m = 1:length(line_depth_2013)
        opt = abs(lat - lat_2013(m)).^2 + abs(lon - lon_2013(m)).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);
        line_depth_2013(m) = depth_2013(row,col);
end

% Calculate the approximate data index of the layer in each trace
surf_ind = round(surf_ind - upper + 1 + line_depth_2013/0.3);

% Convert coordinates to polar stereographic for grid bin finding
[x_2013, y_2013] = ll2psn(lat_2013, lon_2013);

% Calculate the average connectivity in the same 1km grid bins as the 2017
% data
depth_hc_2013 = NaN*ones(1, size(at_axis,1));
mean_lat = zeros(size(at_axis,1),1);
mean_lon = zeros(size(at_axis,1),1);
for m = 1:length(at_axis) - 1
    % Find the 2013 observations that fall in this grid bin
    [~, start13] = min(abs(x_2013 - at_axis(m,1)).^2 + abs(y_2013 - at_axis(m,2)).^2);
    [~, stop13] = min(abs(x_2013 - at_axis(m+1,1)).^2 + abs(y_2013 - at_axis(m+1,2)).^2);
    if stop13 < start13
        a = stop13;
        stop13 = start13;
        start13 = a;
    end
    % Find the center coordinate of the grid bin
    mean_lat(m) = mean(lat_2013(start13:stop13));
    mean_lon(m) = mean(lon_2013(start13:stop13));
    if stop13 > start13 + 10 
        ref = zeros(1, stop13-start13+1);
        meas = zeros(1, stop13-start13+1);
        for k = start13:stop13
            % Total score if there was a perfectly connected layer at that
            % depth
            ref(k-start13+1) = max(perfect_score_2013(surf_ind(k)-8:surf_ind(k)+8,k));
            % Measured value normalized by the perfect score
            meas(k-start13+1) = max(hc_total_2013(surf_ind(k)-8:surf_ind(k)+8,k));
        end
         % Average normalized connectivity score in the grid bin
         depth_hc_2013(m) = nanmean(meas./ref);
    end
end

% New along-track axis
[at_avg, ~, ~, ~] = geodetic_to_along_track(mean_lat, mean_lon, []);

% Final connectivity metric for 2013 data
connectivity_2013 = depth_hc_2013;

%%%%%%%%%%% Display the first overlapping segment %%%%%%%%%%%%%%%%%%%%%%%%%

% Split out the first segment from the combined data
stop = size(at_axis1,1)-1;

[at, ~, ~, ~] = geodetic_to_along_track(mean_lat(1:stop), mean_lon(1:stop), []);
[x1, y1] = ll2psn(mean_lat(1:stop), mean_lon(1:stop));

% Load the MAR-derived total melt production data and find results along
% this transect
load('./RefData/TotalMelt_2013-2017.mat');
[x,y] = ll2psn(lat,lon);
melt = zeros(1,length(x1));
for k = 1:length(melt)
    opt = abs(x - x1(k)).^2 + abs(y - y1(k)).^2;
    [value, index] = min(opt(:));
    [row, col] = ind2sub(size(opt), index);
    melt(k) = melt_holder(row,col);
end

% Save the analysis results used for the figure
seg_lat = mean_lat(1:stop);
seg_lon = mean_lon(1:stop);
density13 = mean_min_density_2013(1:stop);
density17 = mean_min_density_2017(1:stop);
connectivity13 = connectivity_2013(1:stop);
connectivity17 = connectivity_2017(1:stop);
total_melt = melt/1000;
along_track = at/1000;

save('./DerivedData/Demo_TemporalAnalysis_Seg1.mat', 'seg_lat', 'seg_lon', 'density13', 'density17', 'connectivity13', 'connectivity17', 'total_melt', 'along_track');

% Plot the final figure
fig = figure;
left_color = (1/256)*[0 0 0];
right_color = (1/256)*[0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

patch([at'/1000 fliplr(at'/1000)], [0.780*ones(size(at')) 0.830*ones(size(at'))], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on;
plot(at/1000, mean_min_density_2013(1:stop), 'LineWidth', 1, 'Color', (1/256)*[0 114 190]);
hold on;
plot(at/1000, mean_min_density_2017(1:stop), 'LineWidth', 1, 'Color', (1/256)*[218 83 25]);
xlabel('Distance Along-Track (km)');
ylabel('Density (g/cm^3)');
ylim([0 0.92]);
xlim([0 456]);
yyaxis right;
patch(vertcat(at/1000, flipud(at)/1000), [zeros(size(melt)) fliplr(melt)/1000], (1/256)*[240 203 48], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(at/1000, melt/1000, 'Color', (1/256)*[240 203 48], 'LineWidth', 1, 'LineStyle', '-');
hold on;
plot(at/1000, ones(size(at)), '-.', 'Linewidth', 1, 'Color', [0.5 0.5 0.5]);
hold on;
plot(at/1000, connectivity_2013(1:stop), '-', 'LineWidth', 1, 'Color', (1/256)*[96 204 252]);
hold on;
plot(at/1000, connectivity_2017(1:stop), '-', 'LineWidth', 1, 'Color', (1/256)*[235 154 117]); 
ylim([0 2]);
ylabel('Connectivity / Melt (mWe)');
set(gca, 'FontSize', 20, 'FontWeight', 'bold');

%%%%%%%%%%% Display the second overlapping segment %%%%%%%%%%%%%%%%%%%%%%%%

% Split out the second segment from the combined data
start = size(at_axis1,1)+1;

[at, ~, ~, ~] = geodetic_to_along_track(mean_lat(start:end-1), mean_lon(start:end-1), []);

% Load the MAR-derived total melt production data and find results along
% this transect
load('TotalMelt_2013-2017.mat');
[x,y] = ll2psn(lat,lon);
melt = zeros(1, length(mean_lat(start:end-1)));
for k = 1:length(melt)
    [x1, y1] = ll2psn(mean_lat(k + start - 1),mean_lon(k + start - 1));
    opt = abs(x - x1).^2 + abs(y - y1).^2;
    [value, index] = min(opt(:));
    [row, col] = ind2sub(size(opt), index);
    melt(k) = melt_holder(row,col);
end

% Save the results used to plot the figure
seg_lat = mean_lat(start:end-1);
seg_lon = mean_lon(start:end-1);
density13 = mean_min_density_2013(start:end-1);
density17 = mean_min_density_2017(start:end-1);
connectivity13 = connectivity_2013(start:end-1);
connectivity17 = connectivity_2017(start:end-1);
total_melt = melt/1000;
along_track = at/1000;

save('./DerivedData/Demo_TemporalAnalysis_Seg2.mat', 'seg_lat', 'seg_lon', 'density13', 'density17', 'connectivity13', 'connectivity17', 'total_melt', 'along_track');

% Plot the final figure
fig = figure;
left_color = [0 0 0]; 
right_color = [0 0 0]; 
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

patch([at'/1000 fliplr(at'/1000)], [0.780*ones(size(at')) 0.830*ones(size(at'))], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on;
plot(at/1000, mean_min_density_2013(start:end-1), 'LineWidth', 1, 'Color', (1/256)*[0 114 190]);
hold on;
plot(at/1000, mean_min_density_2017(start:end-1), 'LineWidth', 1, 'Color', (1/256)*[218 83 25]);
xlabel('Distance Along-Track (km)');
ylabel('Density (g/cm^3)');
ylim([0 0.92]);
xlim([0 192]);
yyaxis right;
patch(vertcat(at/1000, flipud(at)/1000), [zeros(size(melt)) fliplr(melt)/1000], (1/256)*[240 203 48], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(at/1000, melt/1000, 'Color', (1/256)*[240 203 48], 'LineWidth', 1, 'LineStyle', '-');
hold on;
plot(at/1000, ones(size(at)), '-', 'Linewidth', 1, 'Color', [0.5 0.5 0.5]);
hold on;
plot(at/1000, connectivity_2013(start:end-1), '-', 'LineWidth', 1, 'Color', (1/256)*[96 204 252]);
hold on;
plot(at/1000, connectivity_2017(start:end-1), '-', 'LineWidth', 1, 'Color', (1/256)*[235 154 117]); 
ylim([0 2]);
ylabel({'Connectivity /'; 'Melt (mWe)'});
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
toc