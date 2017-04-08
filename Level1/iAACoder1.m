function [ x ] = iAACoder1( AACSeq1, fNameOut )

y = zeros((length(AACSeq1)+1) * 1024, 2);

for i = 1:length(AACSeq1)
    if strcmp(AACSeq1(i).frameType, 'ESH')
        DCT_frame = zeros(128, 16);
        DCT_frame(:, 1:2:end) = AACSeq1(i).chl.frameF;
        DCT_frame(:, 2:2:end) = AACSeq1(i).chr.frameF;
        
        %temp = iFilterbank(DCT_frame, AACSeq1(i).frameType, AACSeq1(i).winType);
        %y(1024*(i-1) + 449 : 1024*(i-1) + 1600, :) = y(1024*(i-1) + 449 : 1024*(i-1) + 1600, :) + temp;
        
    else
        
        DCT_frame = [AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF];
        %temp = iFilterbank(DCT_frame, AACSeq1(i).frameType, AACSeq1(i).winType);
        %y((i-1)*1024 + 1 : (i+1)*1024, :) = y((i-1)*1024 + 1 : (i+1)*1024, :) + temp;
        
    end
    
    temp = iFilterbank(DCT_frame, AACSeq1(i).frameType, AACSeq1(i).winType);
    y((i-1)*1024 + 1 : (i+1)*1024, :) = y((i-1)*1024 + 1 : (i+1)*1024, :) + temp;
    
end

%We didn't encode the last frame of zeros so we only take out the last 1024
%zeros
y = y(1025:end-1024, :);

if (nargout == 1)
    x = y;
end

audiowrite(fNameOut, y, 48000);

end

