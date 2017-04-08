function [ S, sfc, G, T ] = AACquantizer( frameF, frameType, SMR )
%frameF: 128x8 if ESH, else 1024x1
%frameType: ESH, OLS, LSS, LPS
%SMR: 42x8 if ESH, else 69x1

load('TableB219.mat');
%Use indexing starting at 1
B219a(:, 2:3) = B219a(:, 2:3) + 1;
B219b(:, 2:3) = B219b(:, 2:3) + 1;
MQ = 8191;

if strcmp(frameType, 'ESH')

    temp = cumsum(frameF.^2);
    d = temp(B219b(:, 3), :);
    P = d - [zeros(1, 8); d(1:end-1, :)];
    
    T = P ./ SMR;
    
    %a is 128x8 matrix if ESH
    a = repmat( 16/3 * log2(max(frameF).^(0.75)/MQ), 128, 1 );
    
    mask1 = ones(128, 8);
    mask2 = ones(128, 8);
    acoustic = ones(42, 8);
    diff = ones(42, 8);
    while any(acoustic & diff)
        %Get scalefactor gains for each band
        scfa = a(B219b(:,2), :);
        
        %diff(i) will be 1 if |a(i)-a(i-1)| and |a(i) - a(i+1)| <= 60
        temp  = (abs(scfa - [scfa(2,:); scfa(1:end-1,:)]) <= 60);
        diff = temp .* [temp(2:end,:); temp(end,:)];
        
        for k = 1:42
            mask1(B219b(k,2):B219b(k,3), :) = repmat(diff(k,:), B219b(k,3)-B219b(k,2)+1, 1);
            mask2(B219b(k,2):B219b(k,3), :) = repmat(acoustic(k,:), B219b(k,3)-B219b(k,2)+1, 1);
        end
        
        %If the band's error power s below the threshold and diff(band) is
        %1, then we increment the scalefactor gain
        a = a + mask1.*mask2;
        
        %Calculate the error power for each band
        %No need for the sign() since X_hat has the same sign as X and we
        %square their difference
        e = (frameF - round( (abs(frameF).*2.^(-a/4)).^(3/4) + 0.4054 ).^(4/3) .* 2.^(a/4)).^2;
        temp = cumsum(e.^2);
        d = temp(B219b(:, 3), :);
        Pe = d - [zeros(1, 8); d(1:end-1, :)];
        
        %Compare each band's error power with acoustic threshold
        acoustic = (Pe < T);
    end
    
    %Quantize the MDCT coeffs and the scalefactors; differences
    S = sign(frameF) .* round( (abs(frameF).*(2.^(-a/4))).^(3/4) + 0.4054 );
    sfc = a(B219b(:,3), :) - [zeros(1, 8); a(B219b(1:end-1,3), :)];
    sfc = sign(sfc) .* round( (abs(sfc).*(2.^(-sfc/4))).^(3/4) + 0.4054);
    G   = sfc(1,:);
else
    
    temp = cumsum(frameF.^2);
    d = temp(B219a(:, 3));
    P = d - [0; d(1:end-1)];
    
    T = P ./ SMR;
    
    a = repmat( 16/3 * log2(max(frameF)^0.75 / MQ), 1024, 1 );
    
    mask1 = ones(1024, 1);
    mask2 = ones(1024, 1);
    acoustic = ones(69, 1);
    diff = ones(69, 1);
    while any(acoustic & diff)
        %Get scalefactors for each band
        scfa = a(B219a(:,2));
        
        %Create arrays indicating if difference between bands is > 60
        temp  = (abs(scfa - [scfa(2); scfa(1:end-1)]) <= 60);
        diff = temp .* [temp(2:end); temp(end)];
        
        for k = 1:69
            mask1(B219a(k,2):B219a(k,3)) = diff(k);
            mask2(B219a(k,2):B219a(k,3)) = acoustic(k);
        end
        
        a = a + (mask1.*mask2);
        
        %Calculate power of error signal
        e = (frameF - round( (abs(frameF).*(2.^(-a/4))).^(3/4) + 0.4054 ).^(4/3) .* (2.^(a/4))).^2;
        temp = cumsum(e.^2);
        d = temp(B219a(:, 3));
        Pe = d - [0; d(1:end-1)];
        
        %Compare each band's error power with acoustic threshold
        acoustic = (Pe < T);
    end
    
    %Quantize mdct coeffs and scalefactors
    S = sign(frameF) .* round( (abs(frameF).*(2.^(-a/4))).^(3/4) + 0.4054 );
    sfc = a(B219a(:,3)) - [0; a(B219a(1:end-1, 3))];
    sfc = sign(sfc) .* round( (abs(sfc).*(2.^(-sfc/4))).^(3/4) + 0.4054);
    G   = sfc(1);
    
end


end

