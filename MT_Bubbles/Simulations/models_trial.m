% models_trial
% Thuraisingham 1
% Anderson 2

clc
clear
close all

f_range = linspace(0.1,300,3000)*1000;
a_range = linspace(8e-4,4e-3,10); % very slow for high number of radiuses
% a_range = 4e-3; %a_range(aa); % bubble radius (m)
TS_thur = 10*log10(bubble_response_model(f_range,a_range, 1));
TS_and = 10*log10(bubble_response_model(f_range,a_range, 2));


% Plot freq x TS
a = a_range(1);
figure;
subplot(211)
hold on
plot(f_range/1000, TS_thur);
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
titlename = "Plot freq x TS. TS for a sphere with a=" + (a*100) + " cm";
title(titlename)
subtitlename = "Thuraisingham";
subtitle(subtitlename)

subplot(212)
hold on
plot(f_range/1000, TS_and);
xlabel('Freq (kHz)');ylabel('TS (dB re 1 m^2)')
subtitlename = "Anderson";
subtitle(subtitlename)
%% Plot ka x TS
c=1500;
ka = 2*pi/c*f_range'*a_range;
kk = 1; % at specific radius
figure;
% subplot(211)
semilogx(ka(:,kk), TS_thur(:,kk))
hold on
semilogx(ka(:,kk), TS_and(:,kk))
legend('Thuraisingham','Anderson' )
xlabel('log(ka)');ylabel('TS (dB re 1 m^2)')
titlename = "Plot ka x TS. TS for a sphere with a=" + (a*100) + " cm";
title(titlename)