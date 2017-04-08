function [ frameFout ] = iTNS( frameFin, frameType, TNScoeffs )


if strcmp(frameType, 'ESH')
    r = zeros(4, 8);
    for i = 1:8
        r(:, i) = roots([1; -TNScoeffs(:,i)]);
    end
    if any(any(abs(r) >= 1))
        error('At least one inverse filter is unstable!')
    end
    
    frameFout = zeros(size(frameFin));
    for i = 1:8
        frameFout(:, i) = filter(1, [1 -TNScoeffs(:,i)'], frameFin(:, i));
    end
else
    r = roots([1 -TNScoeffs']);
    if any(abs(r) >= 1)
        error('Inverse filter unstable!')
    end
    
    frameFout = filter(1, [1 -TNScoeffs'], frameFin);
end


end

