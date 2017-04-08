function [ SNR, bitrate, compression ] = demoAAC3( fNameIn, fNameOut, fnameAACoded )

[yin, Fs]  = audioread(fNameIn);

fprintf('=======================\nLevel 3 \n=======================\n');
tic;
AACSeq3 = AACoder3(fNameIn, fnameAACoded);
T1 = toc;
fprintf('Encoding time: %.5f seconds\n', T1);

tic;
yout = iAACoder3(AACSeq3, fNameOut);
T2 = toc;
fprintf('Decoding time: %.5f seconds\n', T2);

yout = yout(1:size(yin, 1), :);

%We use the bitstream as the size of the encoded signal
for k = 1:length(AACSeq3)
    bits = bits + numel([AACSeq3(k).chl.stream AACSeq3(k).chr.stream]);
end

uncofile = dir(fNameIn);
compression = uncofile.bytes*8/bits;
fprintf('Uncompressed audio: %.4f MB\n', uncofile.bytes/2^20);
fprintf('Compressed struct: %.4f KB\n', bits/2^13);
fprintf('Compression ratio: %.4f % (x%.4f)', 100*bits/(uncofile.bytes*8), ...
                                             compression);

                                         
SNR = 10*log10( sum(yin.^2) ./ sum((yin - yout).^2) );
fprintf('Channel 1 SNR: %.5f dB\n', SNR(1));
fprintf('Channel 2 SNR: %.5f dB\n', SNR(2));

bitrate = bits * Fs / size(yin, 1);



end

