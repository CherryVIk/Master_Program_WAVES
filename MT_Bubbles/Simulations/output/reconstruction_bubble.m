clc
clear
close all

load("bubble+signal.mat"); %upload rx with noise, tx from .mat file
rx_nonoise = load("single-bubble-noise.mat").rx; %rx of the single bubble without noise

fs = 192000;
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);

nRxSeqLength = size(rx,1);
NFFT = 32768;
NBins = 16385; 

% Received signal
Rx_noise = fft(rx(end-nSig+1:end, 1), NFFT);
Rx_noise = Rx_noise(1:NBins, :);

Rx_nonoise = fft(rx_nonoise(end-nSig+1:end, 1), NFFT);
Rx_nonoise = Rx_nonoise(1:NBins, :);

% Original signal
% tx signal to freq. domain
Tx = fft(tx, NFFT);
Tx = Tx(1:NBins, :);

% H_hat = R / T
H_hat = Rx_nonoise(:,1) ./ Tx;

tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);

figure;
subplot(211)
plot(tSim(end-nSig+1:end), rx_nonoise(end-nSig+1:end, 1));
grid on;
title('Received signal, time');
subplot(212)
plot(t, tx(:, 1));
grid on;
title('Emitted signal, time');


f = linspace(0, fs/2, NBins);
figure;
subplot(311)
plot(f, Rx_nonoise(:, 1));
grid on;
title('Received non-noisy signal, freq');
subplot(312)
plot(f, Tx(:, 1));
grid on;
title('Emitted signal, freq');
subplot(313)
plot(f, H_hat(:, 1));
grid on;
title('Bubble signal, freq');

