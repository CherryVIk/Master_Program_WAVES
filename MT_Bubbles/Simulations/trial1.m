clear all;
close all;
%% SONAR Parameters
% Speed of sound underwater
cWater = 1500.0;
% Number of transmitting projectors
NTx = 1; % SIMO: Single-input
centerTx = [0 0 0];
% Number of receiving hydrophones
NRx = 32; % SIMO: Multi-output
centerRx = [0 0 0];
% Beamforming angles
angles = -90:2:90;
% Number of beams
NBeams = length(angles);

filename = 'MyAnimation.gif';
moveBubble_step = 10;
bRandomSeed = 42;%42
rng(bRandomSeed)
for time_step = 1:3
%% Signal Parameters
% Sample frequency
fs = 192000;
% Signal bandwidth
fB = 20000;
% Center frequency
fC = 75000;
% Min & max frequencies
fMin = fC - fB/2;
fMax = fC + fB/2;
% Ping duration
tPing = 1.0; % in seconds
nPing = tPing * fs; % in samples
% Signal duration
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);
% Signal types
eSignalTypes.CW = 'CW';
eSignalTypes.blNoise = 'blNoise';
eSignalTypes.HFM = 'HFM';
eSignalType = eSignalTypes.blNoise;


%% Adaption of SONAR-system wrt. signal design
% Set array element distance to half of wavelength of max freq.
dRx = cWater / (fMax) / 2; % Lambda half distance of array elements
dTx = dRx;

%% FFT-parameters
nNextPow2 = nextpow2(nSig*2); % find nearest x so 2^x = nSig
NFFT = 2^nNextPow2; % FFT-length as multiple of 2
% NFFT = nSig;
NBins = NFFT / 2 + 1; % FFT-bins of pos. frequencies
bins = 0:NBins-1; % Freq. bin support vector
fBinRes= fs / NFFT;
nfMin = floor(fMin/fBinRes);
nfMax = ceil(fMax/fBinRes);
f = linspace(-fs/2, fs/2, NFFT);%NFFT
% f = linspace(0, fs/2, NBins);%NFFT

%% Generate transmit sequence
% Bandpass filter design
Fstop1 = fMin-100;       % First Stopband Frequency
Fpass1 = fMin;       % First Passband Frequency
Fpass2 = fMax;       % Second Passband Frequency
Fstop2 = fMax+100;       % Second Stopband Frequency
Astop1 = 100;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
Astop2 = 100;          % Second Stopband Attenuation (dB)
match  = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, fs);
Hd = design(h, 'cheby2');



if strcmp(eSignalType, eSignalTypes.blNoise)
    % Generate Gaussian white noise
    tx = randn(nSig, NTx);
    % an amplitude for that noise is 10% of the noise-free signal at every element.
     % noise_add = 
   
    tx = filter(Hd, tx);
    % tx = tx + noise_add;
    %tx = filtfilt(Hd.sosMatrix, Hd.ScaleValues, tx);
    % Transform time to freq. domain signal
    Tx = fft(tx, NFFT);%NFFT
    % Only save positive freq.
    Tx = Tx(1:NBins, :);
% The following has a commented out ideal (but impractical) bandpass filter
    % Bandpass
    %Tx(1:nfMin, :) = 0;
    %Tx(nfMax:end, :) = 0;
    % Freq. to time domain
    %tx = ifft(Tx, NFFT, 'symmetric');
    %tx = tx(1:nSig);
% End of ideal bandpass example
end

%% Plot transmit sequence
% figure(10);
% subplot(211);
% plot(t, tx(:,1));
% title("Time domain signal");
% grid on;
% subplot(212);
% logTx = 20*log10(abs(Tx(:,1))./max(abs(Tx(:,1))));
% plot(f(end/2:end), logTx);
% title("Log. frequency spectrum ");
% grid on;

%% Bubble environment settings
%  constrains a, b
a=30;
b=35;
Nbubbles=10;
if time_step == 1
    posTar = randi([a,b],Nbubbles,3);
else
    posTarNew = randi([a,b],Nbubbles,3);
    posTar(:,3) = posTar(:,3) + moveBubble_step*ones(NTargets, 1);
    posTar = [posTar; posTarNew];
end
NTargets = size(posTar, 1);
bDirectSound = 0;

x = posTar(:,1);
y = posTar(:,2);
z = posTar(:,3);

figure
plot3(x,y,z, '-o')
zlim([30 60])
grid on
hold on
pause(.5)
% hold on
end