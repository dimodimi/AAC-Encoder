function [ frameT ] = iFilterbank( frameF, frameType, winType )
%frameF: MDCT coefficients 1024x2 if winType = OLS | LSS | LPS
%or 128x16 if winType = ESH

frame = imdct4(frameF);

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
    
    frameT = [win_long win_long] .* frame;
    
elseif (strcmp(frameType, 'LSS'))
    
    lss_win = [win_long(1:1024); ones(448, 1); win_short(129:end); zeros(448, 1)];
    frameT = [lss_win lss_win] .* frame;
    
elseif (strcmp(frameType, 'LPS'))
    
    lps_win = [zeros(448, 1); win_short(1:128); ones(448, 1); win_long(1025:end)];
    frameT = [lps_win lps_win] .* frame;
    
elseif (strcmp(frameType, 'ESH'))
    win_frame = repmat(win_short, 1, 16) .* frame;
    
    frameT = zeros(2048, 2);
    
    for i = 1:8
        frameT(448 + (i-1)*128 + 1:448 + (i+1)*128, :) = ...
            frameT(448 + (i-1)*128 + 1:448 + (i+1)*128, :) + win_frame(:, 2*i-1:2*i);
    end
    
else
    error('Frame type not supported. Accepted values: OLS, LSS, LPS, ESH.')
end


end

