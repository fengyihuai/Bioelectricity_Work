clc; clear all; close all;

for fi = 1:length(frequencies)
    % create wavelet and get its FFT
    fft_wavelet = fft(exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*(wavelet_cycles / (2*pi*frequencies(fi)))^2))/frequencies(fi), n_convolution);
    
    % convolution for all three sites (save only power)
    convolution_result_fft = ifft(fft_wavelet.*fft_data_seed, n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_seed = abs(reshape(convolution_result_fft, EEG.pnts, EEG.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_trgt, n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_trgt = abs(reshape(convolution_result_fft, EEG.pnts, EEG.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_ctrl, n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_ctrl = abs(reshape(convolution_result_fft, EEG.pnts, EEG.trials)).^2;
    
    % downsample and rank transform all data
    conv_result_seed = tiedrank(conv_result_seed(times2saveidx, :)')';
    conv_result_trgt = tiedrank(conv_result_trgt(times2saveidx, :)')';
    conv_result_ctrl = tiedrank(conv_result_ctrl(times2saveidx, :)')';
    
    for ti = 1:length(times2save)
        
        % compute bivariate correlations
        r_st = 1-6*sum((conv_result_seed(ti, :) - conv_result_trgt(ti, :)).^2)/(EEG.trials*(EEG.trials^2-1));
        r_sc = 1-6*sum((conv_result_seed(ti, :)-conv_result_ctrl(ti, :)).^2)/(EEG.trials*(EEG.trials^2-1));
        r_tc = 1-6*sum((conv_result_ctrl(ti, :)-conv_result_trgt(ti, :)).^2)/(EEG.trials*(EEG.trials^2-1));
        
        % bivariate correlation for comparison
        tf_corrdata(fi, ti ,1) = r_st;
        
        % compute partial correalation and store in results matrix
        tf_corrdata(fi, ti, 2) = (r_st-r_sc*r_tc) / ( sqrt(1-r_sc^2)*sqrt(1-r_tc^2) );
    end
end

for i=1:2
    subplot(1, 2, i)
    contourf(times2save, frequencies, squeeze(tf_corrdata(:, :, i)), 40, 'linecolor', 'none')
    set(gca, 'clim', clim, 'xlim', [-200 800], 'yscale', 'log', 'ytick', logspace(log10(frequencies(1)), log10(frequencies(end)), 6), 'yticklabel', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 6)*10)/10)
    axis square
    if i==1
        title(['Correlation between ' seed_chan ' and ' target_chan ])
    else
        title([ 'Partial correlation between ' seed_chan ' and ' target_chan ])
    end
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
end

%% Figure 27.8 

% Re-run the code for the previous figure but comment out the
% following line towards the top;
% times2save = EEG.times; % uncomment this line for figure 27.8
% Then run this section of code.

ds_timesidx = dsearchn(EEG.times', (-200:50:800)'); % ds = downsampled
[~, lofreq] = min(abs(frequencies-4.7));
[~, hifreq] = min(abs(frequencies-32));

figure

subplot(221)
contourf(times2save, frequencies, squeeze(tf_corrdata(:, :, 2)), 40, 'linecolor', 'none')
hold on
plot(get(gca, 'xlim'), frequencies([lofreq lofreq]), 'k--')
plot(get(gca, 'xlim'), frequencies([hifreq hifreq]), 'k--')
set(gca, 'clim', clim, 'xlim', [-200 800], 'yscale', 'log', 'ytick', logspace(log10(frequencies(1)), log10(frequencies(end)), 6), 'yticklabel', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 6)*10)/10)
title('Original (256Hz)')

subplot(222)
contourf(times2save(ds_timesidx), frequencies, squeeze(tf_corrdata(:, ds_timesidx, 2)), 40, 'linecolor', 'none')
hold on
plot(get(gca, 'xlim'), frequencies([lofreq lofreq]), 'k--')
plot(get(gca, 'xlim'), frequencies([hifreq hifreq]), 'k--')
set(gca, 'clim', clim, 'xlim', [-200 800], 'yscale', 'log', 'ytick', logspace(log10(frequencies(1)), log10(frequencies(end)), 6), 'yticklabel', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 6)*10)/10)
title('Down-sampled (20Hz)')

subplot(223)
plot(EEG.times, squeeze(tf_corrdata(lofreq, :, 2)))
hold on
plot(EEG.times(ds_timesidx), squeeze(tf_corrdata(lofreq, ds_timesidx, 2)), 'ro-', 'markerface', 'w')
title('Effect of downsampling on low-frequency activity')
set(gca, 'xlim', [-200 800], 'ylim', [.25 .65])

subplot(224)
plot(EEG.times, squeeze(tf_corrdata(hifreq, :, 2)))
hold on
plot(EEG.times(ds_timesidx), squeeze(tf_corrdata(hifreq, ds_timesidx, 2)), 'ro-', 'markerface', 'w')
title('Effect of downsampling on high-frequency activity')
set(gca, 'xlim', [-200 800], 'ylim', [-.1 .6])
legend({'Original (256 Hz)'; 'Down-sampled (20 Hz)'})

%% Figure 27.9

% note: this cell take a while to run, particularly on slow computers!

n = 1000;
ncorrs = 100000;

t = [0 0 0];
for i=1:ncorrs
    
    % create random variables
    a = rand(2, n);
    
    tic
    % Matlab corr function
    c2 = corr(a', 'type', 's');
    t(1) = t(1) + toc;
    
    
    tic
    % self-written Spearman correlation (must first rank-transform)
    a1 = tiedrank(a')'; % tiedrank accepts matrix input, but make sure the matrix is in the correct orientation!!
    c1 = 1-6*sum((a1(1, :)-a1(2, :)).^2)/(n*(n^2-1));
    t(2) = t(2)+toc;
    
    
    tic
    % ordinary least squares
    % Note: Uncommenting the following line will normalize the data to give
    % you a correlation coefficient, If you don't need the correlation
    % coefficient (and instead can use unstandardized regression
    % coefficients), leave this line commented out for a ten-fold increase
    % in speed. 
    % a = bsxfun(@rdivide, bsxfun(@minus, a, mean(a, 2)), std(a, [], 2));
    c3 = (a(1, :)*a(1, :)')\a(1, :)*a(2, :)';
    t(3) = t(3) + toc;
end

figure
bar(t)
set(gca, 'xticklabel', {'corr function'; 'manual'; 'ols'}, 'xlim', [.5 3.5])
ylabel([ 'Time for ' num2str(ncorrs) ' iterations (s)' ])
