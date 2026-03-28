function filter_norm_audio_fcn(file)

    [wave, fs] = audioread(file);

    %% high-pass filter 
    cutoff = 400; 
    order = 4;
    [b, a] = butter(order, cutoff/(fs/2), 'high');
    wave_filt = filtfilt(b, a, wave);

    %% norm amplitude
    wave_norm = (wave_filt ./ max(abs(wave_filt)))*0.8;
    psychwavwrite(wave_norm, fs, file);

    fprintf('audio filtered (%d high pass)and normalized: %s\n', cutoff, file)

end
