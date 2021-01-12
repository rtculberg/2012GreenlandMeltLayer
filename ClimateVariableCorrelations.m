% Author: Riley Culberg
% Date: 12/7/2020
%
% This script calculates the 2D histograms and distance correlation
% coefficients between layer prominence and various climate variables. This
% script was used to produce Figure 3, Supplementary Figure 3, and Supplementary 
% Table 1.
%
% Outputs:
% Plots and distance correlation coefficients in Figure 3, Supplementary
% Figure 3, and Supplementary Table 1. 
%   x - 1xN vector of doubles, 2D histogram bin edges for climate variables
%   bin_edges - 1xQ vector of doubles, 2D histogram bin edges for layer
%       prominence
%   hist_2d - QxN matrix of doubles, probability of a given layer
%       prominence occuring for a given climate condition
%
% Data dependencies:
% 2012 MARv3.5.2 outputs - /RefData/MAR_DataExtract.mat
% Long term climate averages from MARv3.10 - /RefData/MAR_20km_MeanVals_1980-2011.mat
% Final results of radar analysis - /DemoData/IceLayerDetections/*.txt
%
% Functional Dependencies:
% /ReferenceFunctions/ll2psn.mat - from Arctic Mapping Tools by Chad
% Green
% (https://www.mathworks.com/matlabcentral/fileexchange/63324-arctic-mapping-tools)
% /ReferenceFunctions/distcorr.mat - by Shen Liu
% (https://www.mathworks.com/matlabcentral/fileexchange/39905-distance-correlation)
% ------------------------------------------------------------------------

clear;

% Set up path for reference funtions
addpath(genpath('./ReferenceFunctions/'));

% Load the MAR model data
load('./RefData/MAR_DataExtract.mat');

% Uncomment below for examples of how to run this with the full MAR NetCDF
% files which were WAY too large to package up with the script demos
% lat = ncread('./RefData/ICE_2012_01_12.nc', 'LAT');
% lon = ncread('./RefData/ICE_2012_01_12.nc', 'LON');
% melt_tmp = ncread('./RefData/ICE_2012_01_12.nc', 'ME');
% temp = ncread('./RefData/ICE_2012_01_12.nc', 'TI1');
% surface_temp_2012 = ncread('./RefData/ICE_2012_01_12.nc', 'ST2');
% surface_temp_2011 = ncread('./RefData/ICE.2011.01-12.nc', 'ST2');
% density = ncread('./RefData/ICE_2012_01_12.nc', 'RO1');
% depth = ncread('./RefData/ICE_2012_01_12.nc', 'OUTLAY');
% accum_tmp = ncread('./RefData/ICE_2012_01_12.nc', 'SF');

averages = load('./RefData/MAR_20km_MeanVals_1980-2011.mat');

% Extract the firn density and temperature profiles
density = squeeze(mean(density,4));  % Mean density profile - 2012
temp = squeeze(temp(:,:,:,150));     % Firn temperature profile on 1 June 2012

% Calculate the climate variables of interest
melt = zeros(size(lat,1), size(lat,2));
melt_var = zeros(size(lat,1), size(lat,2));
accum = zeros(size(lat,1), size(lat,2));
summer_temp = zeros(size(lat,1), size(lat,2));
winter_temp = zeros(size(lat,1), size(lat,2));
melt_days = zeros(size(lat,1), size(lat,2));
for k = 1:size(lat,1)
    for m = 1:size(lat,2)
        % Total surface melt production - 2012
        tmp = squeeze(melt_tmp(k,m,1,:));
        tmp2 = squeeze(accum_tmp(k,m,:));
        index = find(tmp >= 0); 
        melt(k,m) = sum(tmp(index));
        % Total annual accumulation - 2012
        accum(k,m) = sum(tmp2);
        % Surface melt rate variability
%%%%%%% !!!!! Comment this line to run the other correlations !!!!!!%%%%%%%
        melt_var(k,m) = std(tmp);
        % Average summer (JJA) surface temperatures in 2012 in degrees K
        summer_temp(k,m) = mean(surface_temp_2012(k,m,1,150:242)) + 273.15;
        % Average winter (DJF) surface temperatures in 2011/2012 in degrees
        % K
        winter_temp(k,m) = 0.5*(mean(surface_temp_2012(k,m,1,1:60)) + mean(surface_temp_2011(k,m,1,334:end))) + 273.15;
        % Number of melt days in 2012
        index = find(tmp > 0);
        melt_days(k,m) = length(index);
    end
end

% Calculate the cold content in the top 20 m of firn on 1 June 2012 and the 
% total latent heat in the total annual surface meltwater production 
% following equation2 in Vandecrux et al (2020) Firn cold content evolution 
% at nine sites on the Greenland ice sheet between 1998 and 2017. 

rho_water = 997;          % density of water at 0C in kg/m^3
heat_capacity = 333.55;   % heat capacity of water in kJ/kg
specific_heat = 2.108;    % specific heat of ice in kJ-kg/K
cold_content = zeros(size(lat,1), size(lat,2));
latent_heat = zeros(size(lat,1), size(lat,2));
for k = 1:size(lat,1)
    for m = 1:size(lat,2)
        
        total_depth = depth(end);
        del_depth = diff(depth);
        
        mass = squeeze(density(k,m,1:end-1)).*del_depth;
        total_mass = sum(mass);
        
        % Depth-weighted average density
        rho_avg = sum(del_depth.*squeeze(density(k,m,1:end-1)))./total_depth;
        % Mass weighted average temperature
        temp_avg = sum(mass.*squeeze(temp(k,m,1:end-1)))./total_mass;
        
        % Cold content in top 20 meters of firn
        cold_content(k,m) = 20*(specific_heat*rho_avg*abs(temp_avg));
        
        % Latent heat
        latent_heat(k,m) = melt(k,m)*((heat_capacity*rho_water)/1000); 
    end
end

% Location of ice layer detection files
base_dir = './DemoData/IceLayerDetections/';

% Flight transects with ice layer detections
% Transects 4-8 are considered the northwest region of Greenland and 9-18
% the southern region
transects = [ 
             "20170328_01" "20170410_01"  "20170413_01" "20170327_04"  "20170412_01"...
             "20170331_01" "20170403_02"...
             "20170414_01" ...
             "20170421_01" "20170422_01" "20170424_01" "20170429_01" ...
             "20170501_02" "20170502_01" "20170505_02" "20170506_01" ...
             "20170508_02" "20170511_01"];

% Convert MAR grid centerpoints to Polar Stereographic coordinates
[x,y] = ll2psn(lat,lon);

variable = [];
layer_prom = [];
for k = 1:18 % edits transects loaded to analyze different regions
    
    % Load ice layer detections
    fprintf('Transect: %s\n', transects(k));
    file = strcat(base_dir, transects(k), '_CleanResults.txt');
    data = readmatrix(file);
    
    % Store the layer prominence metric
    layer_prom = [layer_prom data(:,3)'];
    
    % For every bin in every transect, find the corresponding MAR data
    temp = zeros(length(data(:,3)),1);
    for m = 1:length(data(:,3))
        
        % Find the MAR grid cell nearest to the radar detection
        [x1, y1] = ll2psn(data(m,1),data(m,2));
        opt = abs(x - x1).^2 + abs(y - y1).^2;
        [value, index] = min(opt(:));
        [row, col] = ind2sub(size(opt), index);

%%%%%% Climate variables - uncomment the run the parameter of interest %%%%

%         2012 melt to accumulation ratio
%         temp(m) = melt(row,col)./accum(row,col);

%         Number of melt days in 2012
%         temp(m) = melt_days(row,col);

%         Ratio of 2012 JJA surface temperatures to 2011/2012 DJF
%         temperatures
%         temp(m) = summer_temp(row,col)./winter_temp(row,col);

%         2012 melt anomalies
%         temp(m) = (melt(row,col) - averages.mean_melt(row,col))./averages.std_melt(row,col);

%         Mean annual accumulation
%         temp(m) = averages.mean_accum(row,col);

%         Mean annual melt to accumulation ratio
%         temp(m) = averages.mean_melt(row,col)./averages.mean_accum(row,col);

%         Ratio of latent heat to firn cold content
%         temp(m) = latent_heat(row,col)./cold_content(row,col);

%       Melt variability (standard deviation of daily melt rates)
        temp(m) = melt_var(row,col);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
   
    % Store the climate data
    variable = vertcat(variable, temp);
end

% Calculate the distance correlation and display results
% correlation = distcorr(variable, layer_prom');
% fprintf('Distance Correlation: %f\n', correlation);

% Plot the 2D histograms showing correlation between layer prominence and
% the climate variable of interest

bin_edges = 0:0.01:1;  % vertical (layer prominence) bin size
x = min(variable):0.03:max(variable);  % horiztonal (climate variable) bin size

% X spacing settings used in manuscript plots 
% 0.01 for 2012 M/A
% 0.03 for melt variability
% 0.005 for cold content
% 10 for mean annual accumulation
% 0.003 for mean annual melt to accumulation ratio 
% 0.15 for 2012 melt anomalies 
% 0.04 for 2012 z-scored anomalies 
% 0.0005 for temperature swing

hist_2d = zeros(length(bin_edges), length(x));
num = 0;
for k = 1:length(x) - 1
    ind = [];
    for m = 1:length(variable)
        if variable(m) >= x(k) && variable(m) < x(k+1)
            ind = [ind m];
        end
    end
    tmp1 = variable(ind);
    tmp2 = layer_prom(ind);
    data_slice = [];
    % Sum up number of observations that fall within each 2D bin
    for m = 1:length(bin_edges)-1
        count = 0;
        for q = 1:length(tmp2)
            if tmp2(q) >= bin_edges(m) && tmp2(q) < bin_edges(m+1)
                count = count + 1;
                data_slice = [data_slice tmp2(q)];
            end
        end
        % Normalize total observations in a given 1D slice
        hist_2d(m,k) = count/length(tmp1);
    end
end

% Plot 2D histogram figure
figure;
imagesc(x, bin_edges, hist_2d);
set(gca, 'YDir', 'normal');
h = colorbar;
caxis([0 0.15]);
cmocean('tempo');
xlabel('2012 Melt Variability');
ylabel('Layer Prominence');
set(get(h,'label'),'string','Probability');
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
xlim([0.2 6]);  
ylim([0.1 1]);
axis square;

save('./DerivedData/Demo_ClimateCorrelations.mat', 'x', 'bin_edges', 'hist_2d');

% Xlim setting used in manuscript plots:
% [0.06 2.5] for 2012 M/A
% [0.2 6] for melt variability
% [0 1] for cold content
% [0 2000] for mean annual accumulation
% [0 1] for mean annual melt to accumulation ratio
% [0 20] for 2012 z-scored anomalies
% [1.07 1.136] for temperature swing
