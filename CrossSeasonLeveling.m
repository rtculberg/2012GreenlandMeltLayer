% Author: Riley Culberg
% Date: 12/7/2020
%
% This scripts uses cleaned surface echo powers to geostatistically
% cross-level a set of flight transects. It gives a set of dc offsets that
% minimize the difference in surface power at crossovers in the mean square
% sense. This is not an absolute calibration, but a relative calibration
% that accounts for system differences between flights. 
% 
% Output:
%   dc_offsets - Mx1 double vector of calibration constants in dB where M
%   is the number of flight transects being leveled and the order of the
%   constants is the same order in which the flight transects occur in the
%   input flight_tracks data structure
%
% Data dependencies:
% Cleaned surface power data - /DemoData/SurfaceData_2017.mat
%
% Function dependencies:
% ./ReferenceFunctions/find_xovers.m
% ./ReferenceFunctions/find_fine_xovers.m
% ./ReferenceFunctions/ll2psn.m
% CVX - this script uses CVX for the least squares leveling which can be 
% downloaded here: http://cvxr.com/cvx/download/
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Reference list for flight_tracks data structure:
% 1 - 20120330_01
% 2 - 20120404_01
% 3 - 20120411_01
% 4 - 20120413_01
% 5 - 20120418_01
% 6 - 20120419_01
% 7 - 20130402_01
% 8 - 20130405_01
% 9 - 20130415_01
% 10 - 20170328_01
% 11 - 20170422_01
% 12 - 20170424_01
% 13 - 20170327_04
% 14 - 20170331_01
% 15 - 20170403_02
% 16 - 20170410_01
% 17 - 20170412_01
% 18 - 20170413_01 
% 19 - 20170414_01
% 20 - 20170421_01
% 21 - 20170429_01
% 22 - 20170501_02
% 23 - 20170502_01
% 24 - 20170505_02
% 25 - 20170506_01
% 26 - 20170508_02
% 27 - 20170511_01

load('./DemoData/Demo_SurfaceData_2017.mat');
flight_tracks = SurfaceData(1:end);
clear SurfaceData;

% Choose the subset of flight tracks to use for leveling (for example, this
% paper levels the northern and southern transects separately to the same
% tie network of transects 1-4 & 6-7 & 9-10)

flight_tracks = {flight_tracks{1:4} flight_tracks{6:7} flight_tracks{9:10} flight_tracks{13:19}};

% Change from geographic to stereographic polar coordinates
for k = 1:length(flight_tracks)
    [x,y] = ll2psn(flight_tracks{k}.lat, flight_tracks{k}.lon);
    flight_tracks{k}.x = x;
    flight_tracks{k}.y = y;
end

% Calculate the crossover tolerance based on the average distance between
% traces across all transects
spacing = 0;
for k = 1:length(flight_tracks)
    spacing = 0.5*(mean(abs(diff(flight_tracks{k}.x))) + mean(abs(diff(flight_tracks{k}.y)))) + spacing;
end

tolerance = spacing/length(flight_tracks);

% To improve efficiency, calculate rough crossovers first to localize the search area 
surveys = 1:1:length(flight_tracks);
num_xover_guess = 6000;     % approximate # of crossovers expected
mean_pow_radius = 50;       % distance over which to average surface power
search_area = 100;          % factor by which to downsample the trace spacing

% Downsample flight tracks
downsampled_flight_tracks = flight_tracks;
for k = 1:length(flight_tracks)
    downsampled_flight_tracks{k}.x = downsampled_flight_tracks{k}.x(1:search_area:end);
    downsampled_flight_tracks{k}.y = downsampled_flight_tracks{k}.y(1:search_area:end);
    downsampled_flight_tracks{k}.surf_geo_pow = downsampled_flight_tracks{k}.surf_geo_pow(1:search_area:end);
end

% Find the rough crossover points
xovers = find_xover(surveys, downsampled_flight_tracks, search_area*tolerance, num_xover_guess, search_area*mean_pow_radius);

% Find the average distance between traces
sampling = 0;
for k = 1:length(flight_tracks)
    sampling = sampling + mean(sqrt(diff(flight_tracks{k}.x).^2 + diff(flight_tracks{k}.y).^2));
end
sampling = sampling/length(flight_tracks);

% Calculate the fresnel radius of the system and set the distance over
% which to average the surface power (here we use 10 times the Fresnel
% radius)
lambda = 3e8/725e6;          % Radar wavelength
h = 500;                     % Approximate flight altitude
rf = sqrt((lambda*h)/2);     % Radius of first fresnel zone
rf = round(10*rf/sampling);  % Averaging distance

% Calculate exact crossover points localized by rough crossover
% approximations
fine_xovers = find_fine_xover(flight_tracks, xovers, search_area, rf);

% Setup data structure for leveling
power = zeros(length(fine_xovers),2);
survey = zeros(length(fine_xovers),2);
for k = 1:size(power, 1)
    survey(k, 1) = fine_xovers{k}.survey1;
    survey(k, 2) = fine_xovers{k}.survey2;
    power(k, 1) = fine_xovers{k}.survey1_pow;
    power(k,2) = fine_xovers{k}.survey2_pow;
end

matches.power = power;
matches.survey = survey;

% Least squares leveling of crossover powers
disp(['Uncorrected RMSD: ' ...
        num2str(norm(matches.power(:,1)-matches.power(:,2)) / ...
                sqrt(size(matches.survey,1)))])
cvx_begin quiet
    variable dc_offset(max(matches.survey(:)),1)
    variable adj_pows1(size(matches.survey,1),1) 
    variable adj_pows2(size(matches.survey,1),1)
    %use huber penalty function, kinked at 3dB
    minimize (norm(adj_pows1 - adj_pows2))
    subject to
        adj_pows1 == matches.power(:,1) + dc_offset(matches.survey(:,1))
        adj_pows2 == matches.power(:,2) + dc_offset(matches.survey(:,2))  
cvx_end
assert(strcmp(cvx_status, 'Solved'))
disp(['Corrected RMSD: ' ...
        num2str(norm(adj_pows1-adj_pows2) / ...
                sqrt(size(matches.survey,1)))])

% Normalize offsets to tie transect that is used for absolute calibration
% to firn cores
dc_offset = dc_offset - dc_offset(1);
            
% Save results
save('./DerivedData/Demo_CalibrationConstants.mat', 'dc_offset');
