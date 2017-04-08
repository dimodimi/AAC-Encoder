function [ AACSeq2 ] = AACoder2( fNameIn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

winType = 'KBD';
[y, Fs] = audioread(fNameIn);
if mod(size(y, 1), 1024)
    y = ([zeros(1024, 2); y; zeros(3072 - mod(size(y, 1), 1024), 2)]);
else
    y = ([zeros(1024, 2); y; zeros(2048, 2)]);
end

%Assume frame -2 is OLS
prevFrameType = 'OLS';

field1 = 'frameType'; 
field2 = 'winType';   value2 = winType;
field3 = 'chl';
field4 = 'chr';

for i = 1:1024:size(y, 1)-3071
    
    value1 = SSC(y(i:i+2047, :), y(i+1024:i+3071, :), prevFrameType);
    prevFrameType = value1;
    
    frameF = filterbank(y(i:i+2047, :), value1, winType);
    
    if (strcmp(value1, 'ESH'))
        [frameTNS, TNScoeffs] = TNS(frameF(:, 1:2:end), value1);
        value3 = struct('frameF', frameTNS, 'TNScoeffs', TNScoeffs);
        
        [frameTNS, TNScoeffs] = TNS(frameF(:, 2:2:end), value1);
        value4 = struct('frameF', frameTNS, 'TNScoeffs', TNScoeffs);
    else
        [frameTNS, TNScoeffs] = TNS(frameF(:, 1), value1);
        value3 = struct('frameF', frameTNS, 'TNScoeffs', TNScoeffs);
        
        [frameTNS, TNScoeffs] = TNS(frameF(:, 2), value1);
        value4 = struct('frameF', frameTNS, 'TNScoeffs', TNScoeffs);
    end
    
    AACSeq2((i-1)/1024 + 1) = struct(field1, value1, field2, value2, field3, value3, field4, value4);
    
end

end

