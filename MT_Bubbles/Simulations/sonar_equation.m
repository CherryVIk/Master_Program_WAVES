% Date: 02.04.2024
% Name: sonar equation implementation
% Wave propagation simulation from sonar to the bubble with the TS calculation
clc; clear; close all

SL = 220; % dB re 1mPa at 1m, source level
AG = 20; % dB, array gain

R = 40; % m, range, radius
TLg = 20*log10(R);% spherical TL from geometrical spreading, one way
alpha_att = 0.5; % dB/km, absorption coefficient
TLl = alpha_att * R; % dB, attenuation ~ TL from dissipation
TL = TLg + TLl; % dB, transmission loss 

r0 = 4e-3; % bubble radius
TS = 10*log10(r0^2/4); % dB, target strength
NL = 60; % dB, noise level
T = 0.1; % s, pulse length
BW = 1/T; % bandwidth of the receiver
NL_tot = NL + 10*log10(BW);

SNR = SL - 2*TL + TS - (NL - AG); % signal to noise ration, active sonar equation 
%% 1. Calculate the transmission loss from SONAR to 1m in front of target
function TL = TLsonar_to1m_tar(R, alpha_att)
    TLg = 20*log10(R-1);% spherical TL from geometrical spreading, one way
    TLl = alpha_att * (R-1); % dB, attenuation ~ TL from dissipation
    TL = TLg + TLl; % dB, transmission loss 
end
%% 2. Calculate TL from 1m in front of target to target
function TL = TLtar_to1m_tar(alpha_att)
    TLg = 20*log10(1);% spherical TL from geometrical spreading, one way
    TLl = alpha_att * (1); % dB, attenuation ~ TL from dissipation
    TL = TLg + TLl; % dB, transmission loss 
end
%% 3. Calculate TS from filtering with frequency response of target
function TS = TStar(sigma_bs)
    TS = 10*log10(sigma_bs); 
end
%% 4. Add transmission loss from target to 1 m in front of target
% - Combine these for total target strength
function TL = TS_total(R, alpha_att, sigma_bs) 
    TL = (TLsonar_to1m_tar + TLtar_to1m_tar) * 2 + 
end
