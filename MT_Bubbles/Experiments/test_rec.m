clear all
close all 
clc

a = audioread('recordings/ping_037_internal_ppfHydIn_16.wav');
X = a(:,1);
%Fs = 5000;%marDeCangas
%Fs = 500000;%whistle
Fs = 192000;
Y = fft(X);
L = length(Y);
Y = abs(Y/L);
% Y = Y(1:L/2+1); 
logY = Y(1:L/2+1);
% logY = 20*log10(abs(Y(:,1))./max(abs(Y(:,1))));

xlen = length(X);
t = (0: xlen-1) / Fs;     
f = Fs*(0:(L/2))/L;
figure
subplot(2,1,1);plot(t,X)
subplot(2,1,2);plot(f,logY) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|Y(f)|')
figure
spectrogram(X,2048,2048/2,L,Fs,'yaxis');
colormap("jet")
 saveas(gcf,"internal_ppfHydIn_15"+".png");

%% Range: [a,b] 
% find at 33-34 s, 47-48, 47.2-47.3
a = 0; b = 1;
% a = 47.2; b = 47.3;
[x, tmin] = min(abs(t-a));
[x, tmax] = min(abs(t-b));
% X_range = X(tmin:tmax);
X_range = X;
Y_range = fft(X_range);
L_y = length(Y_range);
Y_range = abs(Y_range/L);
Y_range = Y_range(1:L_y/2+1); 
f_range = Fs*(0:(L_y/2))/L_y;

figure;
% subplot(2,1,1);plot(t(tmin:tmax),X_range)
subplot(2,1,1);plot(t,X_range)
subplot(2,1,2);plot(f_range,Y_range)
title('Single-Sided Amplitude Spectrum of X(t), ranged')
xlabel('f (Hz)')
ylabel('|Y(f)|')
% subplot(2,1,2);plot(f_range,20*log10(abs(Y_range(:,1))./max(abs(Y_range(:,1))))) 