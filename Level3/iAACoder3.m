function [ x ] = iAACoder3( AACSeq3, fNameOut )
%AACSeq3: array of structs of encoded signal
%fNameOut: path to save decoded signal

y = zeros((length(AACSeq3)+1)*1024, 2);

for i = 1:length(AACSeq3)
    
    if strcmp(AACSeq3(i).frameType, 'ESH')
       
        S = [decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook) ...
             decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook)];
         
        sfc = [decodeHuff(AACSeq3(i).chl.sfc, AACSeq3(i).chl.codebook, 12) ...
               decodeHuff(AACSeq3(i).chr.sfc, AACSeq3(i).chr.codebook, 12)];
           
        frameTNS = zeros(128, 16);
        frameTNS(:, 1:2:end) = iAACquantizer(S(:,1), sfc(:,1), ...
                                AACSeq3(i).chl.G, AACSeq3(i).chl.frameType);
        frameTNS(:, 2:2:end) = iAACquantizer(S(:,2), sfc(:,2), ...
                                AACSeq3(i).chr.G, AACSeq3(i).chr.frameType);
                            
        DCTframe = zeros(128, 16);
        DCTframe(:, 1:2:end) = iTNS(frameTNS(:, 1:2:end), AACSeq3(i).frameType, ...
                                    AACSeq3(i).chl.TNScoeffs);
        DCTframe(:, 2:2:end) = iTNS(frameTNS(:, 2:2:end), AACSeq3(i).frameType, ...
                                    AACSeq3(i).chr.TNScoeffs);
        
    else
        
        S = [decodeHuff(AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook) ...
             decodeHuff(AACSeq3(i).chr.stream, AACSeq3(i).chl.codebook)];
         
        sfc = [decodeHuff(AACSeq3(i).chl.sfc, AACSeq3(i).chl.codebook, 12) ...
               decodeHuff(AACSeq3(i).chr.sfc, AACSeq3(i).chr.codebook, 12)];
           
        frameTNS = [iAACquantizer(S(:,1), sfc(:,1), AACSeq3(i).chl.G, AACSeq3(i).chl.frameType) ...
                    iAACquantizer(S(:,2), sfc(:,2), AACSeq3(i).chr.G, AACSeq3(i).chr.frameType)];
                
        DCTframe = [iTNS(frameTNS(:,1), AACSeq2(i).frameType, AACSeq2(i).chl.TNScoeffs) ... 
            iTNS(frameTNS(:,2), AACSeq2(i).frameType, AACSeq2(i).chr.TNScoeffs)];
        
    end
    
    temp = iFilterbank(DCTframe, AACSeq3(i).frameType, AACSeq3(i).winType);
    y((i-1)*1024 + 1 : (i+1)*1024, :) = y((i-1)*1024 + 1 : (i+1)*1024, :) + temp;
    
end

y = y(1025:end-1024, :);

if (nargout == 1)
    x = y;
end

audiowrite(fNameOut, y, 48000);

end

