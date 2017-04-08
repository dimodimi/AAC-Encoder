function [ SMR ] = psycho( frameT, frameType, frameTprev1, frameTprev2 )
%frameT : mono sound - 2048 x 1 matrix
%frmaType : OLS, LSS, ESH, LPS
%frameTprev1 : 256 x 1 if frameType else 2048 x 1
%frameTprev2 : same

load('TableB219.mat');
load('spreads.mat');

%Use indexing starting at 1
B219a(:, 2:3) = B219a(:, 2:3) + 1;
B219b(:, 2:3) = B219b(:, 2:3) + 1;

if strcmp(frameType, 'ESH')
    
    %s: 256 x 10 matrix - 2 previous + 8 subframes
    s = [frameTprev2 frameTprev1 buffer(frameT(449:end-448, 1), 256, 128)];
    s = s(:, [1:2 4:10]);
    
    sw = repmat(0.5 - 0.5*cos(pi*(0:255 + 0.5)/256), 1, 10) .* s;
    
    sfr = fft(sw);
    r = abs(sfr(1:128, :));
    f = angle(sfr(1:128, :));
    
    rpred = 2 * r(:, 2:9) - r(:, 1:8);
    fpred = 2 * f(:, 2:9) - f(:, 1:8);
    r = r(:, 3:10);
    f = f(:, 3:10);
        
    c = sqrt( (r.*cos(f) - rpred.*cos(fpred)).^2 + (r.*sin(f) - rpred.*sin(fpred)).^2 ) ...
        ./ (r + abs(rpred));

    temp = cumsum([r c] .* [r r.*r]);
    d = temp(B219b(:, 3), :);
    d2 = d - [zeros(1, 16); d(1:end-1, :)];
    e  = d2(:,1:8);
    c2 = d2(:,9:16);
    
    ecb = reshape( sum( kron(e, ones(1, 42)) .* repmat(spread_short, 1, 8) ), 42, 8 );
    ct  = reshape( sum( kron(c2, ones(1, 42)) .* repmat(spread_short, 1, 8) ), 42, 8 );
    
    cb = ct ./ ecb;
    en = ecb ./ repmat(sum(spread_short)', 1, 8);
    
    tb = -0.299 - 0.43*log(cb);
    
    SNR = 6 * tb + 18 * (1 - tb);
    
    bc = 10 .^ (-SNR/10);
    
    nb = en .* bc;
    
    npart = max( nb, repmat(eps()*128*(10.^(B219b(:, end)/10)), 1, 8) );
    
    SMR = e ./ npart;
    
else
    
    sw = repmat(0.5 - 0.5*cos(pi*(0:2047 + 0.5)/2048)', 1, 3) .* [frameTprev2 frameTprev1 frameT];
    
    temp = fft(sw);
    r = abs(temp(1:1024, :));
    f = angle(temp(1:1024, :));
    
    rpred = 2*r(:, 2) - r(:, 1);
    fpred = 2*f(:, 2) - f(:, 1);
    r = r(:, 3);
    f = f(:, 3);
    
    c = sqrt((r.*cos(f) - rpred.*cos(fpred)).^2 + (r.*sin(f) - rpred.*sin(fpred)).^2) ...
        ./ (r + abs(rpred));
    
    temp = cumsum([r c] .* [r r.*r]);
    d = temp(B219a(:, 3), :);
    d2 = d - [0 0; d(1:end-1, :)];
    e = d2(:,1);
    c2 = d2(:,2);
    
    ecb = sum(repmat(e, 1, 69) .* spread_long)';
    ct  = sum(repmat(c2, 1, 69) .* spread_long)';
    
    cb = ct ./ ecb;
    en = ecb ./ (sum(spread_long)');
    
    tb = -0.299 - 0.43*log(cb);
    
    SNR = 6 * tb + 18 * (1 - tb);
    
    bc = 10 .^ (-SNR/10);
    
    nb = en .* bc;
    
    npart = max( nb, eps*1024*(10.^(B219a(:, end)/10)) );
    
    SMR = e ./ npart;
    
end

end

