function [layers, layer_depth] = LayerSimulator_UWB(fc, BW, tau, n, n_depth, sampling, bin, width, ice_density)
    c = 3e8;
    
    f = (fc - BW/2):10e6:(fc + BW/2);

    core_diff = diff(n_depth);
    R = zeros(length(f), length(n_depth));
    distance = zeros(1, length(n_depth));
    index = zeros(1, length(n_depth));

    R_index = 1;
    ind1 = 1;  
    depth = n_depth(ind1);
    saved_depth = depth;
    total_t = 0;
    delta_t = [];
    tracker = 1;
    
    for m = 1:length(f)
        lambda = c/f(m);
        R_index = 1;
        ind1 = 1;
        depth = n_depth(ind1);
        tracker = 1;
        while depth < n_depth(end)
            distance(1, R_index) = depth;
            
            del_t = 0;
            ind2 = ind1;
            while del_t < (1/(2*BW))
                if ind2 > length(n) - 1
                    break;
                end
                del_t = del_t + core_diff(ind2)*(n(ind2)/c);
                ind2 = ind2 + 1;
            end
            
            if ind2 > length(n) - 1
                ind2 = length(n) - 1;
            end
            
            if tracker == bin
                middle = round(0.5*(ind2 - ind1) + ind1);
                space = floor(0.5*width);
                n(middle - space:middle + space) = 1 + 0.845*ice_density;
            end
            
            depth = depth + sum(core_diff(ind1:ind2));
            
            if m == 1
                total_t = total_t + del_t;
                delta_t = [delta_t del_t];
                saved_depth = [saved_depth depth];
            end
            
            index(1, R_index) = ind1;
            
            tmp = TransferMatrix(lambda, 0, n(ind1:ind2+1), core_diff(ind1:ind2-2));
            R(m, R_index) = tmp(1);
            R_index = R_index + 1;
            ind1 = ind2;
            tracker = tracker + 1;
        end
        
    end
    
    R = R(:,1:R_index-1);
    distance = distance(:,1:R_index-1);
    distance = distance + (c/(4*BW*1.6));
    index = index(1,1:R_index-1);
    
    Fs = 3e9;           % Sampling frequency (Hz)
    Fc = 0;
    slope = BW/tau;     % Chirp slope
    
    chirp = SingleSideBand(slope, tau, Fs, Fc);
    time_window = tukeywin(length(chirp),0.2);
    chirp = chirp.*time_window';
    
    freq_axis = linspace(0, Fs/2, ceil(length(chirp)/2));
    [~, ind] = min(abs(freq_axis - BW));
    freq_win = hann(ind)';
    freq_win = [freq_win zeros(1,length(chirp) - ind)];
    
    ref = fft(chirp, length(chirp));
    ipr = fftshift(ifft(conj(ref).*fft(chirp).*freq_win));
    
    reflectivity = zeros(1, size(R,2));
    for m = 1:size(R,2)
        freq_reflec = interp1(f - (fc - BW/2), R(:,m), freq_axis);
        
        freq_reflec(isnan(freq_reflec)) = mean(abs(R(:,m)));
        freq_reflec = [freq_reflec mean(abs(R(:,m)))*ones(1, length(ref) - length(freq_reflec))];
        
        reflec_spectrum = ref.*freq_reflec;
        pc = ifft(fftshift(conj(ref).*reflec_spectrum.*freq_win));
        reflectivity(m) = max(abs(pc).^2)./max(abs(ipr).^2);
    end
        
    full_t = 0:(1/Fs):(2*total_t);
    delta_reflec = zeros(size(full_t));
    
    next = 0;
    count = 1;
    while next < full_t(end) && count <= length(reflectivity)
        [~, ind] = min(abs(full_t - next));
        delta_reflec(ind) = reflectivity(count);
        next = next + 2*delta_t(count);
        count = count + 1;
    end
    
    result = conv((abs(ipr).^2)./max(abs(ipr).^2), delta_reflec);
    result2 = result((floor(length(ipr)/2)-1):(end - floor(length(ipr)/2) + 1));
    dist = (0:1:length(result2)-1)*(1/Fs)*(c/(2*1.6));
    
    x = 0:sampling:dist(end);
    resamp = interp1(dist, result2, x);
    resamp(1) = max(result2(1:50));
    
    layers = resamp;
    layer_depth = x;
    
end