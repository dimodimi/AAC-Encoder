function [ AACSeq3 ] = AACoder3( fNameIn, fnameAACoded )
%Create array of structure for each frame
%For each channel we store the frame type, window type, TNS coeffs,
%acoustic threshold, global gain, quantized scalefactors, bitsream
%of quantized MDCT coeeficients and the huffman codebook used for 
%the MDCT coeffs quantization.

huffLUT = loadLUT();
winType = 'KBD';
[y, Fs] = audioread(fNameIn);
if mod(size(y, 1), 1024)
    y = ([zeros(3072, 2); y; zeros(3072 - mod(size(y, 1), 1024), 2)]);
else
    y = ([zeros(3072, 2); y; zeros(2048, 2)]);
end

prevFrameType = 'OLS';

field1 = 'frameType'; 
field2 = 'winType';   value2 = winType;
field3 = 'chl';
field4 = 'chr';

for i = 2049:1024:size(y, 1)-3071
    
    value1 = SSC(y(i:i+2047, :), y(i+1024:i+3071, :), prevFrameType);
    prevFrameType = value1;
    
    frameF = filterbank(y(i:i+2047, :), value1, winType);
    
    if strcmp(value1, 'ESH')
        [frameTNS, TNScoeffs] = TNS(frameF(:, 1:2:end), value1);
        SMR = psycho(y(i:i+2047, 1), value1, y(i-832:i-577, 1), y(i-704:i-449, 1));
        [S, sfc, G, T] = AACquantizer(frameTNS(:, 1:2:end), value1, SMR);
        [huffmdct, huffcodebook] = encodeHuff(S, huffLUT);
        [huffsfc, ~]   = encodeHuff(sfc, huffLUT, 12);
        
        value3 = struct('TNScoeffs', TNScoeffs, 'T', T, 'G', G, ...
            'sfc', huffsfc, 'stream', huffmdct, 'codebook', huffcodebook);
        
        [frameTNS, TNScoeffs] = TNS(frameF(:, 2:2:end), value1);
        SMR = psycho(y(i:i+2047, 2), value1, y(i-832:i-577, 2), y(i-704:i-449, 2));
        [S, sfc, G, T] = AACquantizer(frameTNS(:, 2:2:end), value1, SMR);
        [huffmdct, huffcodebook] = encodeHuff(S, huffLUT);
        [huffsfc, ~]  = encodeHuff(sfc, huffLUT, 12);
        
        value4 = struct('TNScoeffs', TNScoeffs, 'T', T, 'G', G, ...
            'sfc', huffsfc, 'stream', huffmdct, 'codebook', huffcodebook);
    else
        
        [frameTNS, TNScoeffs] = TNS(frameF(:, 1), value1);
        SMR = psycho(y(i:i+2047, 1), value1, y(i-1024:i+1023, 1), y(i-2048:i-1, 1));
        [S, sfc, G, T] = AACquantizer(frameTNS(:, 1), value1, SMR);
        [huffmdct, huffcodebook] = encodeHuff(S, huffLUT);
        [huffsfc, ~] = encodeHuff(sfc, huffLUT, 12);
        
        value3 = struct('TNScoeffs', TNScoeffs, 'T', T, 'G', G, ...
            'sfc', huffsfc, 'stream', huffmdct, 'codebook', huffcodebook);
        
        [frameTNS, TNScoeffs] = TNS(frameF(:, 2), value1);
        SMR = psycho(y(i:i+2047, 2), value1, y(i-1024:i+1023, 2), y(i-2048:i-1, 2));
        [S, sfc, G, T] = AACquantizer(frameTNS(:, 2), value1, SMR);
        [huffmdct, huffcodebook] = encodeHuff(S, huffLUT);
        [huffsfc, ~] = encodeHuff(sfc, huffLUT, 12);
        
        value4 = struct('TNScoeffs', TNScoeffs, 'T', T, 'G', G, ...
            'sfc', huffsfc, 'stream', huffmdct, 'codebook', huffcodebook);
    end
    
    AACSeq3((i-1)/1024 - 1) = struct(field1, value1, field2, value2, field3, value3, field4, value4);
end

save(fnameAACoded, AACSeq3);

end

