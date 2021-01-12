% Author: Riley Culberg
% Date: 12/7/2020
%
% This function calculates the radar response to a firn density profile and
% allows the user to specify the the locations where ice layers should be
% added.
%
% Inputs:
%   fc - center frequency of radar system in Hz (scalar double)
%   FW - bandwidth of radar system in Hz (scalar double)
%   tau - radar pulse length in s (scalar double)
%   n - permittivity profile of the subsurface (M x 1 vector double)
%   n_depth - depth axis for the permittivity profile (M x 1 vector double)
%   sampling - radar sample spacing in m (scalar double)
%   bin - range bin in which to add an ice layer (scalar double)
%   width - ice layer thickness in # of permittivity samples (scalar double)
%   ice_density - density of ice layer in g/cm^3 (scalar double)
%
% Outputs:
%   layers - reflectivity in linear power units for each range bin (Qx1
%   vector double)
%   layer_depth - depth in meters of each reflectivity (Qx1 vector double)
% 
% Function dependencies:
% ./ReferenceFunctions/TransferMatrix.m
% ./ReferenceFunctions/SingleSideband.m
% ------------------------------------------------------------------------

function [layers, layer_depth] = LayerSimulator_UWB(fc, BW, tau, n, n_depth, sampling, bin, width, ice_density)
    
    c = 3e8;  % speed of light in a vacuum
    
    f = (fc - BW/2):10e6:(fc + BW/2);  % frequency sampling

    core_diff = diff(n_depth);  % Sample spacing in permittivity profile
    
    % Variables to save results
    R = zeros(length(f), length(n_depth));
    distance = zeros(1, length(n_depth));
    index = zeros(1, length(n_depth));
    
    % Start from the beginning of the permittivity profile
    R_index = 1;
    ind1 = 1;  
    depth = n_depth(ind1);
    saved_depth = depth;
    total_t = 0;
    delta_t = [];
    tracker = 1;
    
    % For each frequency step
    for m = 1:length(f)
        lambda = c/f(m);
        R_index = 1;
        ind1 = 1;
        depth = n_depth(ind1);
        tracker = 1;
        % Over the top 15 meters of the firn
        while depth < 15
            distance(1, R_index) = depth;
            
            % Find all of the permittivity sample which fall within the
            % current range bin
            del_t = 0;
            ind2 = ind1;
            while del_t < (1/(2*BW))
                if ind2 > length(n) - 1
                    break;
                end
                del_t = del_t + core_diff(ind2)*(n(ind2)/c);
                ind2 = ind2 + 1;
            end
            
            % Make sure we don't got beyond the end of the permittivity
            % profile
            if ind2 > length(n) - 1
                ind2 = length(n) - 1;
            end
            
            % If we're in the range bin where the ice layer should be
            % added, add it at the center of the bing with the desired
            % width and density
            if tracker == bin
                middle = round(0.5*(ind2 - ind1) + ind1);
                space = floor(0.5*width);
                n(middle - space:middle + space) = 1 + 0.845*ice_density;
            end
            
            % Increment the depth to the end of this range bin
            depth = depth + sum(core_diff(ind1:ind2));
            
            % Increment the two-way travel time axis to end of range bin
            if m == 1
                total_t = total_t + del_t;
                delta_t = [delta_t del_t];
                saved_depth = [saved_depth depth];
            end
            
            index(1, R_index) = ind1;
            
            % Calculate the effective reflection coefficient of the stack
            % of layers in the range bin
            tmp = TransferMatrix(lambda, 0, n(ind1:ind2+1), core_diff(ind1:ind2-2));
            R(m, R_index) = tmp(1);
            % Move to next range bin
            R_index = R_index + 1;
            ind1 = ind2;
            tracker = tracker + 1;
        end
        
    end
    
    % Clip data to region where we actually calculated stuff
    R = R(:,1:R_index-1);
    distance = distance(:,1:R_index-1);
    distance = distance + (c/(4*BW*1.6));
    index = index(1,1:R_index-1);
    
    Fs = 3e9;           % Sampling frequency (Hz)
    Fc = 0;             % Baseband center frequency
    slope = BW/tau;     % Chirp slope
    
    % Model the transmitted signal from the radar
    chirp = SingleSideBand(slope, tau, Fs, Fc);
    time_window = tukeywin(length(chirp),0.2);
    chirp = chirp.*time_window';
    
    % Center the frequency windowing function on the transmit signal spectrum in
    % the frequency domain
    freq_axis = linspace(0, Fs/2, ceil(length(chirp)/2));
    [~, ind] = min(abs(freq_axis - BW));
    freq_win = hann(ind)';
    freq_win = [freq_win zeros(1,length(chirp) - ind)];
    
    % Calculate the impulse response
    ref = fft(chirp, length(chirp));
    ipr = fftshift(ifft(conj(ref).*fft(chirp).*freq_win));
    
    % Calculate the true reflectivity by weighting the frequency domain 
    % transmit spectrum by the reflectivity as a function of frequency,
    % pulse compress the signal, and take the peak reflectivity 
    reflectivity = zeros(1, size(R,2));
    for m = 1:size(R,2)
        freq_reflec = interp1(f - (fc - BW/2), R(:,m), freq_axis);
        
        freq_reflec(isnan(freq_reflec)) = mean(abs(R(:,m)));
        freq_reflec = [freq_reflec mean(abs(R(:,m)))*ones(1, length(ref) - length(freq_reflec))];
        
        reflec_spectrum = ref.*freq_reflec;
        pc = ifft(fftshift(conj(ref).*reflec_spectrum.*freq_win));
        reflectivity(m) = max(abs(pc).^2)./max(abs(ipr).^2);
    end
    
    % Full two-way travel time axis
    full_t = 0:(1/Fs):(2*total_t);
    delta_reflec = zeros(size(full_t));
    
    % Set up a comb function in depth/time space with the layer
    % reflectivities
    next = 0;
    count = 1;
    while next < full_t(end) && count <= length(reflectivity)
        [~, ind] = min(abs(full_t - next));
        delta_reflec(ind) = reflectivity(count);
        next = next + 2*delta_t(count);
        count = count + 1;
    end
    
    % Convolve with the impulse response of the system
    result = conv((abs(ipr).^2)./max(abs(ipr).^2), delta_reflec);
    result2 = result((floor(length(ipr)/2)-1):(end - floor(length(ipr)/2) + 1));
    dist = (0:1:length(result2)-1)*(1/Fs)*(c/(2*1.6));
    
    % Resample to the desired depth sampling 
    x = 0:sampling:dist(end);
    resamp = interp1(dist, result2, x);
    resamp(1) = max(result2(1:50));
    
    % Output data
    layers = resamp;
    layer_depth = x;
    
end