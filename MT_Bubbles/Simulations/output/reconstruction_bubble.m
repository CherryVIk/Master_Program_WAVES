clc
clear
close all

% load("bubble+signal.mat"); %upload rx with noise, tx from .mat file
load("tx_signal.mat"); %upload tx from .mat file
% trans_sig = load("rx-no-bubble-resp.mat").rx; %rx of the single target without bubble response
bubble_sig = load("rx-with-bubble-resp.mat").rx; %rx of the single target with bubble response

fs = 192000;
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);

nRxSeqLength = size(bubble_sig,1);
NFFT = 32768;
NBins = 16385; 

% Transmitted signal
X = fft(tx, NFFT);
X = X(1:NBins, :);

% Received signal
Y = fft(bubble_sig, NFFT);
Y = Y(1:NBins, :);

% H_hat = R / T
H_hat = Y(:,1) ./ X(:,1);
H_hat = abs(H_hat);

tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);

figure;
subplot(311)
plot(t, tx);
grid on;
title('Target signal, time');
subplot(312)
plot(tSim, bubble_sig);
grid on;
title('Bubble signal, time');


f = linspace(0, fs/2, NBins);
figure;
subplot(311)
% logRx_noB = X(:,1);
logRx_noB = 20*log10(abs(X(:,1))./max(abs(X(:,1))));
plot(f, logRx_noB(:, 1));
grid on;
title('Target signal, freq');
subplot(312)
% logRx_withB = Y(:,1);
logRx_withB = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));
plot(f, logRx_withB(:, 1));
grid on;
hold on
title('Target-bubble signal, freq');
subplot(313)
% logH_hat = H_hat(:,1);
logH_hat = 20*log10(abs(H_hat(:,1))./max(abs(H_hat(:,1))));
plot(f, logH_hat(:, 1));
title('Extracted bubble signal, freq');

