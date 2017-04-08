function [ SNR ] = demoAAC1( fNameIn, fNameOut )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[yin, Fs]  = audioread(fNameIn);

fprintf('=======================\nLevel 1 \n=======================\n');
tic;
AACSeq = AACoder1(fNameIn);
T1 = toc;
fprintf('Encoding time: %.5f seconds\n', T1);
tic;
yout = iAACoder1(AACSeq, fNameOut);
T2 = toc;
fprintf('Decoding time: %.5f seconds\n', T2);
yout = yout(1:size(yin, 1), :);

SNR = 10*log10( sum(yin.^2) ./ sum((yin - yout).^2) );

fprintf('Channel 1 SNR: %.5f dB\n', SNR(1));
fprintf('Channel 2 SNR: %.5f dB\n', SNR(2));

end

