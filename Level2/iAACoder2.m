function [ x ] = iAACoder2( AACSeq2, fNameOut )

y = zeros((length(AACSeq2)+1) * 1024, 2);

for i = 1:length(AACSeq2)
    if strcmp(AACSeq2(i).frameType, 'ESH')
%         DCT_frame_TNS = zeros(128, 16);
%         DCT_frame_TNS(:, 1:2:end) = AACSeq2(i).chl.frameF;
%         DCT_frame_TNS(:, 2:2:end) = AACSeq2(i).chr.frameF;
        
        DCT_frame = zeros(128, 16);
        DCT_frame(:, 1:2:end) = iTNS(AACSeq2(i).chl.frameF, AACSeq2(i).frameType, AACSeq2(i).chl.TNScoeffs);
        DCT_frame(:, 2:2:end) = iTNS(AACSeq2(i).chr.frameF, AACSeq2(i).frameType, AACSeq2(i).chr.TNScoeffs);
        
%         temp = iFilterbank(DCT_frame, AACSeq2(i).frameType, AACSeq2(i).winType);
%         y(:, 1024*(i-1) + 449 : 1024*(i-1) + 1600) = y(:, 1024*(i-1) + 449 : 1024*(i-1) + 1600) + temp;
%         
    else
        
%         DCT_frame_TNS = [AACSeq2(i).chl.frameF AACSeq2(i).chr.frameF];
        DCT_frame = [iTNS(AACSeq2(i).chl.frameF, AACSeq2(i).frameType, AACSeq2(i).chl.TNScoeffs) ... 
            iTNS(AACSeq2(i).chr.frameF, AACSeq2(i).frameType, AACSeq2(i).chr.TNScoeffs)];
        
%         temp = iFilterbank(DCT_frame, AACSeq2(i).frameType, AACSeq2(i).winType);
%         y(:, (i-1)*1024 + 1 : (i+1)*1024) = y(:, (i-1)*1024 + 1 : (i+1)*1024) + temp;
%         
    end
    
    temp = iFilterbank(DCT_frame, AACSeq2(i).frameType, AACSeq2(i).winType);
    y((i-1)*1024 + 1 : (i+1)*1024, :) = y((i-1)*1024 + 1 : (i+1)*1024, :) + temp;
end

y = y(1025:end-1024, :);

if (nargout == 1)
    x = y;
end

audiowrite(fNameOut, y, 48000);

end

