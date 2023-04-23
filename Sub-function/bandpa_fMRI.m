function [bandpasSig, hilbertSig, ampSig, phaseSig] = bandpa_fMRI(sigValid,fs,fLow,fHigh) ;
% bandpass the fMRI data
% output: bandpassSig, hilbertSig, ampSig, phaseSig

% % Design Butterworth band-pass filter
% h = fdesign.bandpass('N,F3dB1,F3dB2',6,fLow,fHigh,fs);
% Hd = design(h, 'butter');
% set(Hd, 'Arithmetic', 'double');
%     
% SOS = Hd.sosMatrix;
% G = Hd.ScaleValues;

[z,p,k]=butter(4,[fLow fHigh]/(fs/2),'bandpass');
% [z,p,k]=butter(4,[fLow]/(fs/2),'low');
[sos,g]=zp2sos(z,p,k);
[b, a]=butter(4,[fLow fHigh]/(fs/2),'bandpass');
% [b, a]=butter(4,[fLow]/(fs/2),'low');

badChannels = find(isnan(sigValid(:,1))==1) ;
bandpasSig = nan(size(sigValid)) ;
hilbertSig = nan(size(sigValid)) ;
phaseSig = nan(size(sigValid)) ;
ampSig = nan(size(sigValid)) ;

% which filtfilt
disp('start fitlering...');
tic
for iNode = setdiff(1:size(sigValid,1),badChannels)
    sigOriDemean = (sigValid(iNode,:)-mean(sigValid(iNode,:))) ;
    % bandpasSig(iNode,:) = filter(b,a,sigValid(iNode,:));
    bandpasSig(iNode,:) = filtfilt(sos,g,sigOriDemean);   % may be corrupted by
    % fieldtrip toolbox
    hilbertSig(iNode,:) = hilbert(bandpasSig(iNode,:)) ;
    phaseSig(iNode,:) = angle(hilbertSig(iNode,:)) ;
    ampSig(iNode,:) = abs(hilbertSig(iNode,:)) ;
end
toc