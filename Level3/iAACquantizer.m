function [ frameF ] = iAACquantizer( S, sfc, G, frameType )
%Returns the MDCT coeffs after applying the inverse quantization

%Get the scalefactor gains from the quantized differences
%These will have some accumulated quantization error
sf = cumsum(sfc);

if strcmp(frameType, 'ESH')
    a = zeros(1024, 8);
    for b = 1:42
        a(B219b(b,2):B219b(b,3), :) = repmat(sf(b, :), B219b(b,3)-B219b(b,2)+1,1);
    end
    
    frameF = sign(S) .* (abs(S) .^ (4/3)) .* (2 .^ (a/4));
else
    a = zeros(1024, 1);
    for b = 1:69
        a(B219a(b,2):B219a(b,3), :) = sf(b);
    end
    
    frameF = sign(S) .* (abs(S) .^ (4/3)) .* (2 .^ (a/4));
end

end

