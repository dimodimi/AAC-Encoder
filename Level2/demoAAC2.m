function [ SNR ] = demoAAC2( fNameIn, fNameOut )

[yin, Fs]  = audioread(fNameIn);

fprintf('=======================\nLevel 2 \n=======================\n');
tic;
AACSeq = AACoder2(fNameIn);
T1 = toc;
fprintf('Encoding time: %.5f seconds\n', T1);

tic;
yout = iAACoder2(AACSeq, fNameOut);
T2 = toc;
fprintf('Decoding time: %.5f seconds\n', T2);

yout = yout(1:size(yin, 1), :);

SNR = 10*log10( sum(yin.^2) ./ sum((yin - yout).^2) );

fprintf('Channel 1 SNR: %.5f dB\n', SNR(1));
fprintf('Channel 2 SNR: %.5f dB\n', SNR(2));

end

