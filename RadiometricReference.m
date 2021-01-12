% Author: Riley Culberg
% Date: 12/7/2020
%
% This script calculates the expected radar reflectivity of an icy layer
% givens its density, thickness, and the density of the background firn.
% This calculation uses the transfer matrix based method described in
% Culberg & Schroeder, 'Strong Potential for the Detection of Ice Layers in
% Greenland's Firn by Airborne Radar Sounding', IGARSS 2020 to calculate
% the effective reflection coefficient of a stack of layers in a range bin.
%
% Outputs:
%   layer_reflectivity - an MxN matrix of doubles representing the mean
%   layer reflectivity in dB. M is the range of layer thickness and N is
%   the range of layer densities. 
%
% Data dependencies:
% Fitting coefficients for the mean density profiles -
%    ./RefData/mean_density_model.mat
% Fitting coefficients for the ARMA model describing the seasonal density
% variability - ./RefData/coefficients.mat
% Range of seasonal variability - ./RefData/variability.mat
%
% Function dependencies:
% ./ReferenceFunctions/LayerSimulater_UWB.m
% ./ReferenceFunctions/TransferMatrix.m
% ------------------------------------------------------------------------

clear;

addpath(genpath('./ReferenceFunctions/'));

% Mean density profile
load('./RefData/mean_density_model.mat');
% Autoregressive random process model coefficients used to stochastically
% simulate the seasonal density variability (trained on high resolution
% core data presented in Horhold, et al 2011 - see Culberg & Schroeder,
% 2020b for details of the stochastic density model)
load('./RefData/coefficients.mat'); 
% Range of variability envelope (see Horhold, 2011)
load('./RefData/variability.mat');

% Additional data needed for ARMA model
coefficients.Constant = 3.0051e-5;
coefficients.Variance = 0.0106;

Fc = 750e6;        % Radar center frequency (Hz)
BW = 300e6;        % Radar bandwidth (Hz)
tau = 2.048e-6;    % Radar pulse length (s)
sampling = 0.3;    % Effective sample spacing (m)

% Range bin - this is the depth at which the ice layer will be simulated as
% so it effectively set the background density. To simulate at a particular
% density, find the depth at which that density is reach in the mean
% density profile, then divide that depth by the sample spacing (see
% above). 
bin = 15;          

del_x = 0.003;     % Permittivity model depth resolution (m)

a = model.a;       % Mean density profile model coefficients
b = model.b;
c = model.c;
d = model.d;

% Generate deterministic components of example density profile
density_depth = 0:del_x:12;
mean_density = a*exp(b*density_depth) + c*exp(d*density_depth);
env_model = fit([density_depth(1) density_depth(end)]', [var1 var2]', 'exp1');
envelope = env_model.a*exp(env_model.b*density_depth);

% Layer thickness and density parameter space
thickness = 0.01:0.005:0.3;
ice_density = 0.6:0.01:0.92;
% Number of random density profile trials to run for each parameter
% combination
num_trials = 10;

% Set up variables to hold results
layer_reflectivity = zeros(length(thickness), length(ice_density));
layer_variability = zeros(length(thickness), length(ice_density));
surface_reflectivity = zeros(length(thickness), length(ice_density));
surface_variability = zeros(length(thickness), length(ice_density));
tic
for m = 1:length(thickness)
    fprintf('Thickness: %f\n', thickness(m));
    for p = 1:length(ice_density)
        
        layer = zeros(1, num_trials);
        surface = zeros(1, num_trials);
        for q = 1:num_trials
            % Generate random density profile
            variability = simulate(coefficients, length(density_depth));
            density = mean_density + (envelope.*variability');
            % Convert density to permittivty using Kovacs relationship
            n = 1 + 0.845*density;
            n = [1 n];
            n_depth = [density_depth density_depth(end) + del_x];
            
            % Simulate the the radar response to the density profile
            [layers, layer_depth] = LayerSimulator_UWB(Fc, BW, tau, n, n_depth, sampling, bin, thickness(m)/del_x, ice_density(p));
            
            % Find the layer reflectivity
            [~, ind] = min(abs(layer_depth - bin*sampling));
            [tmp, ~] = max(layers(ind-3:ind+3));
            
            layer(q) = tmp;
            surface(q) = layers(1);
        end
        
        % Save the mean reflectivity of all trials
        layer_reflectivity(m,p) = mean(layer);
        
    end
    fprintf('Percent: %f\n', 100*m/length(thickness));
end
toc

% Save results
save('./DerivedData/LayerReflectRef.mat', 'layer_reflectivity');
