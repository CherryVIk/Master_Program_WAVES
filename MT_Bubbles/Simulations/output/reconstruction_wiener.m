clc
clear
close all

load("bubble+signal.mat"); %upload rx, tx from .mat file
rx_noise = load("single-bubble+noise.mat").rx; %rx
rx_nonoise = load("single-bubble-noise.mat").rx; %rx

fs = 192000;
tSig = 0.05; % in seconds
nSig = tSig * fs; % in samples
t = linspace(0, tSig, nSig);

nRxSeqLength = size(rx,1);
NFFT = 32768;
NBins = 16385; 

% Received signal
Rx_noise = fft(rx_noise, NFFT);
Rx_noise = Rx_noise(1:NBins, :);

Rx_nonoise = fft(rx_nonoise, NFFT);
Rx_nonoise = Rx_nonoise(1:NBins, :);

% Original signal
% tx signal to freq. domain
Tx = fft(tx, NFFT);
Tx = Tx(1:NBins, :);

tSim = linspace(0, nRxSeqLength/fs, nRxSeqLength);
figure;
subplot(211)
plot(tSim, rx_noise(:, 1));
grid on;
title('Received signal, time');
subplot(212)
plot(t, tx(:, 1));
grid on;
title('Emitted signal, time');

f = linspace(0, fs/2, NBins);
figure;
subplot(311)
plot(f, Rx_noise(:, 1));
grid on;
title('Received noisy signal, freq');
subplot(312)
plot(f, Rx_nonoise(:, 1));
grid on;
title('Received non-noisy signal, freq');
subplot(313)
plot(f, Tx(:, 1));
grid on;
title('Emitted signal, freq');


%% Goal: obtain the frequency response of the single bubble
% y = rx(end-nSig+1:end, 1); 
% x = tx; 
y = rx_nonoise(:,1); x = rx_noise(:,1); 
y = y(:); % reference signal
x = x(:); % signal with additive Gaussian noise
N = 7000; % filter order
[xest,b,MSE] = wienerFilt(x,y,N);
%% plot results

figure
subplot(311)
plot(tSim(N+1:end),x(N+1:end),'r')
hold on
plot(tSim(N+1:end),y(N+1:end),'k')
ylim([-1e-5,1e-5]);
title('Wiener filtering example')
legend('noisy signal','reference')
subplot(312)
plot(tSim(N+1:end),xest,'k')
ylim([-1e-5,1e-5]);
legend('estimated signal')
subplot(313)
plot(tSim(N+1:end),(x(N+1:end) - xest),'k')
legend('residue signal')
xlabel('time (s)')
ylim([-1e-5,1e-5]);

%% Functions
% https://de.mathworks.com/matlabcentral/fileexchange/71440-signal-separation-with-wiener-filtering
   
function [xest,B,MSE] = wienerFilt(x,y,N)
    X = 1/N .* fft(x(1:N));
    Y = 1/N .* fft(y(1:N));
    X = X(:);
    Y = Y(:);
    Rxx = N .* real(ifft(X .* conj(X))); % Autocorrelation function
    Rxy = N .* real(ifft(X .* conj(Y))); % Crosscorrelation function
    Rxx = toeplitz(Rxx);
    Rxy = Rxy';
    B = Rxy / Rxx; B = B(:); % Wiener-Hopf eq. B = inv(Rxx) Rxy
    xest = fftfilt(B,x);
    xest = xest(N+1:end); % cut first N samples due to distorsion during filtering operation
    MSE = mean(y(N+1:end) - xest) .^2; % mean squared error
end