clear all
close all 
clc

a = audioread('recordings/internal_ppfHydIn_15.wav');
X = a(:,1);
%Fs = 5000;%marDeCangas
%Fs = 500000;%whistle
Fs = 192000;
Y = fft(X);
L = length(Y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure
subplot(2,1,1);plot(X)
subplot(2,1,2);plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
figure
spectrogram(X,1280,1200,1280,Fs,'yaxis');
 saveas(gcf,"whistle16x_Fs_500000"+".png");