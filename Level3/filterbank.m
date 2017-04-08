function [ frameF ] = filterbank( frameT, frameType, winType )
%frameT: 2048x2 matrix
%frameType: OLS, LSS, ESH or LPS
%winType: KBD (Kaiser) or SIN (sinusoidal)

if (strcmp(winType, 'KBD'))
    
    klong  = kaiser(1025, 6);
    kshort = kaiser(129, 4);
    
    acc_long  = cumsum(klong(1:1024));
    acc_short = cumsum(kshort(1:128));
    
    win_long  = sqrt( [acc_long; acc_long(end:-1:1)] / sum(klong) );
    win_short = sqrt( [acc_short; acc_short(end:-1:1)] / sum(kshort) );
    
elseif (strcmp(winType, 'SIN'))
    win_long  = sin(pi*((0:2047) + 0.5)/2048)';
    win_short = sin(pi*((0:255) + 0.5)/256)';
else
    error('Window type not supported. Accepted values: KBD or SIN')
end

if (strcmp(frameType, 'OLS'))
    
    win_frame = [win_long win_long] .* frameT;
    
elseif (strcmp(frameType, 'LSS'))
    
    lss_win = [win_long(1:1024); ones(448, 1); win_short(129:end); zeros(448, 1)];
    win_frame = [lss_win lss_win] .* frameT;
    
elseif (strcmp(frameType, 'LPS'))
    
    lps_win = [zeros(448, 1); win_short(1:128); ones(448, 1); win_long(1025:end)];
    win_frame = [lps_win lps_win] .* frameT;
    
elseif (strcmp(frameType, 'ESH'))
    
    ch1_frames = buffer(frameT(449:end-448, 1), 256, 128);
    ch2_frames = buffer(frameT(449:end-448, 2), 256, 128);
    frame = zeros(256, 16);
    frame(:, 1:2:end) = ch1_frames(:, 2:end);
    frame(:, 2:2:end) = ch2_frames(:, 2:end);
    
    win_frame = repmat(win_short, 1, 16) .* frame;
    
else
    error('Frame type not supported. Accepted values: OLS, LSS, LPS, ESH.')
end

frameF = mdct4(win_frame);


end

