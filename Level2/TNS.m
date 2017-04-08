function [ frameFout, TNScoeffs ] = TNS( frameFin, frameType )
%FrameFin: MDCT coeeficients - 128x8 for ESH else 1024x1 (1 channel)
%P: 
%TNScoefffs: Quantized TNS coeffs - 4x8 for ESH else 4x1

load('TableB219.mat');

%Array indexing starts with 1
b_long  = B219a(:, 2:3)+1;
b_short = B219b(:, 2:3)+1;

if strcmp(frameType, 'ESH')
    %If frameType = ESH then P = 128x8 matrix
    S = zeros(128, 8);
    for i = 1:size(b_short, 1)
        S(b_short(i,1):b_short(i,2), :) = bsxfun(@times, ...
            sqrt(sum(frameFin(b_short(i,1):b_short(i,2), :).^2)), ones(b_short(i,2)-b_short(i,1)+1, 8));
    end
    
    for k = 127:-1:1
        S(k, :) = (S(k,:) + S(k+1,:))/2;
    end
    for k = 2:128
        S(k,:) = (S(k,:) + S(k-1,:))/2;
    end

%     S = cumsum( cumsum(S ./ repmat(2.^[1:127 127]', 1, 8), 'reverse') ...
%         .* repmat( (2.^([0:127]- [127 127:-1:1]))', 1, 8)) .* repmat( (2.^[127:-1:0])', 1, 8 );
    
else
    S = zeros(1024, 1);
    for i = 1:size(b_long, 1)
        S(b_long(i,1):b_long(i,2), :) = ...
            sqrt(sum(frameFin(b_long(i,1):b_long(i,2), :).^2)) * ones(b_long(i,2)-b_long(i,1)+1, 1);
    end
        
    for k = 1023:-1:1
        S(k) = (S(k) + S(k+1))/2;
    end
    for k = 2:1024
        S(k) = (S(k) + S(k-1))/2;
    end
    
%     S = cumsum( cumsum(S ./ 2.^[1:1023 1023]', 'reverse') ...
%         .* (2.^([0:1023]-[1023 1023:-1:1]))' ) .* (2.^[1023:-1:0])';
end

X = frameFin ./ S;

if strcmp(frameType, 'ESH')
    R = zeros(32);
    p = zeros(32, 1);
    for i = 1:8
        [r, lags] = xcorr(X(:,i), 4);
        r = r(lags >= 0);
        
        R((i-1)*4 + 1 : 4*i, (i-1)*4 + 1 : 4*i) = toeplitz(r(1:end-1));
        p(4*(i-1)+1:4*i) = r(2:end);
    end
    TNScoeffs = reshape(R \ p, 4, 8);
else
    [r, lags] = xcorr(X, 4);
    r = r(lags >= 0);
    R = toeplitz(r(1:end-1));
    
    TNScoeffs = R \ r(2:end);
end

%Quantize
TNScoeffs(TNScoeffs <= -0.7) = -0.75;
TNScoeffs(TNScoeffs > 0.7) = 0.75;
TNScoeffs((TNScoeffs > -0.7) & (TNScoeffs <= 0.7)) = ...
    floor(10*TNScoeffs((TNScoeffs > -0.7) & (TNScoeffs <= 0.7)))/10 + 0.05;

%Filter
if strcmp(frameType, 'ESH')
    frameFout = zeros(128, 8);
    for i = 1:8
        frameFout(:, i) = filter([1 -TNScoeffs(:,i)'], 1, frameFin(:, i));
    end
else
    frameFout = filter([1 -TNScoeffs'], 1, frameFin);
end

end

