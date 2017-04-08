function [ AACSeq1 ] = AACoder1( fNameIn )
%fName1: Wav filename - 2 channels - 48 kHz sampling rate

winType = 'KBD';
[y, Fs] = audioread(fNameIn);
if mod(size(y, 1), 1024)
    y = [zeros(1024, 2); y; zeros(3072 - mod(size(y, 1), 1024), 2)];
else
    y = [zeros(1024, 2); y; zeros(2048, 2)];
end

%Assume frame -2 is OLS
prevFrameType = 'OLS';

field1 = 'frameType'; 

field2 = 'winType';   value2 = winType;
field3 = 'chl';       
field4 = 'chr';

%We don't encode the last frame of 2048 zeros.
for i = 1:1024:size(y, 1)-3071
    
    value1 = SSC(y(i:i+2047, :), y(i+1024:i+3071, :), prevFrameType);
    prevFrameType = value1;
    
    frameF = filterbank(y(i:i+2047, :), value1, winType);
    
    if (strcmp(value1, 'ESH'))
        value3 = struct('frameF', frameF(:, 1:2:end));
        value4 = struct('frameF', frameF(:, 2:2:end));
    else
        value3 = struct('frameF', frameF(:, 1));
        value4 = struct('frameF', frameF(:, 2));
    end
    
    AACSeq1((i-1)/1024 + 1) = struct(field1, value1, field2, value2, field3, value3, field4, value4);
    
end

end

