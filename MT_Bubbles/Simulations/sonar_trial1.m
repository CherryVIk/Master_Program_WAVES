%%% Author: VB
%%% Start Date: 1.02.2024
%%% update

clear
close all

c0 = 1500;
fs = 192000;
fb = 2e4;
fc = 4e4;

nTr = 1;
nRec = 10;

lambda = c0/(fc+fb/2);
dR = lambda/2;
angles = -90:5:90;
nAngles = length(angles);

T = 0.1;
dt = 1/fs;
t = 0:dt:T;
nT = length(t);

% filter
Fst1 = 1.99e4; %First Stopband Frequency
Fp1 = 2e4;%First Passband Frequency
Fp2 = 5e4;
Fst2 = 5.01e4;
Ast1 = 100;  %attenuation, dB
Ap = 1; % ripple, fluctuations in the passband, dB
Ast2 = 100;
cheb_filter = fdesign.bandpass(Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2, fs);
Hd = design(cheb_filter, 'ellip');

% Noise
noise = randn(nT, nTr);
noise = filter(Hd, noise);

nF = 2^nextpow2(nT);
fBinRes = fs / nF;
f = linspace(0,fs/2,nF/2+1); %to get real-valued signal
fft_noise = fft(noise,nF);
fft_noise = fft_noise(1:nF/2+1);

% Plotting
figure;
subplot(211)
plot(t,noise);title('Noise, time')
log_noise = 20*log10(abs(fft_noise)./max(abs(fft_noise)));
subplot(212)
plot(f,log_noise);title('Noise, log freq')

% Targets
posTar = [0 40; 30 40];
nTar = size(posTar,1);

center = [0 0];
centrRec = (nRec+1)/2;
centrTr = (nTr+1)/2;
angleRec = [cosd(0) sind(0)];
angleTr = [cosd(0) sind(0)];
posRec = (((1:nRec) - centrRec) * dR).'*angleRec + center;
posTr = (((1:nTr) - centrRec) * dR).'*angleTr + center;

% Calculate delays
distProp = zeros(nTr,nRec,nTar);
for iTr = 1:nTr
    for iRec = 1:nRec
        for iTar = 1:nTar
            dTrToPos = norm(posTar(iTar,:)-posTr(iTr,:));
            dRecToPos = norm(posTar(iTar,:)-posRec(iTr,:));
            distProp(iTr,iRec, iTar) = dTrToPos + dRecToPos;
        end
    end
end
timeProp = distProp./c0; %m / m/s = s
timeProp = timeProp.*fs; % s * 1/s = 1

% Calculate received signals
nSeq = nT + ceil(max(timeProp(:))); % ???
recSignals = zeros(nSeq, nRec);
for iTr = 1:nTr
    for iRec = 1:nRec
        for iTar = 1:nTar
            indS = floor(timeProp(iTr,iRec, iTar))+1; % ???
            indE = floor(indS + length(noise(:,iTr)))-1; % ???
            recSignals(indS:indE,iRec) = recSignals(indS:indE,iRec) + noise(:,iTr);
%             HERE for bubble response
        end
    end
end
propTxSample = (0:nSeq)/fs;
propDxSample = propTxSample * c0;

%% Correlation
% we calculate correclation for the centered receiver/transmitter
[corr,lags] = xcorr(recSignals(:,round(nRec/2)), noise(:,round(nTr/2)));
corr = abs(corr(nSeq:end))./max(abs(corr(nSeq:end))); % Only consider positive correlation lags
tlagsMeters = lags(nSeq:end)./fs.*c0;
% find peaks
[pks, locs] = findpeaks(corr, 'MinPeakHeight', 0.8, 'MinPeakDistance',6);
% locsLeast = min(squeeze(timeProp(round(nTr/2),:, round(nRec/2))));
timeSim = linspace(0,nSeq/fs,nSeq);

figure;
subplot(211)
plot(timeSim, recSignals(:,1));title('Received signal')
subplot(212);
plot(tlagsMeters, corr);
hold on;
plot(tlagsMeters(locs), pks, 'rx');
title('Crosscorrelation: Transmit- & receive signal');

%% Beamforming: AMV
nFft = 2^nextpow2(2*nSeq);
nBins = nFft / 2 + 1;
bins = 0:nBins - 1;
f = bins * fBinRes;
tTau = (sin(angles) .* posRec(:,1))/c0;
F = repmat(f, [nAngles, 1, nRec]);
F = permute(F, [1 3 2]);
tauXf = bsxfun(@times,tTau',F);
AMV = exp(1j*2*pi*tauXf);
AMV = permute(AMV, [2 3 1]);
%% Matched filtering 
RecSignals = fft(recSignals,nFft);
RecSignals = RecSignals(1:nBins,:);
Noise = fft(noise,nFft);
Noise = Noise(1:nBins,:);
Mf = conj(Noise) .* RecSignals; % Matched filtering: noise seq in received seq
MfBf = squeeze(sum(Mf.' .* AMV, 1)).'; % Delay-and-sum
mfbf = ifft(MfBf, nFft, 2, 'symmetric'); %time domain
mfbf = abs(mfbf(:,1:nSeq));
nSize = 1000;

%% Downscale
downscale = nSeq / nSize;
mfbfS = zeros(nAngles, nSize);
for ii = nSize:-1:1
    indS = round((ii-1)*downscale + 1);
    indE = round(ii*downscale);
    mfbfS(:,ii)=max(mfbf(:,indS:indE),[],2);
end

%% Plot PPI (plan position indicator)
figure;
ppi = 20*log10(abs(mfbfS)+eps);
% ppi(ppi<-30)=-30;
dMax = nSeq/fs*c0/2;
theta = deg2rad(angles + 90);
r = 0:dMax/size(ppi,2):dMax-1/size(ppi,1);
[Theta, Rr] = meshgrid(theta, r);
[A,B] = pol2cart(Theta, Rr);
surf(A,B,ppi.');
colormap('jet');
view(0,90);
xlabel('x [m]');
ylabel('y [m]');
daspect([1 1 1]);
axis tight
shading interp;
colorbar;
ax = gca;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf, 'units', 'pixels', 'position', [100 40 1500 900]);
hold on;
hold off;