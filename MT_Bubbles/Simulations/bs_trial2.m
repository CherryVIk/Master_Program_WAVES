clear 
close all

% f = 10e4; 
% a = 50e-6;
f_range = linspace(10e3,1e6,1000); % echosounder freq (Hz=1/s)
a_range = linspace(1e-3,20e-3,1000);  % bubble radius (m)
%a_range = linspace(6e-6,2e-4,1000);  % bubble radius (m)
%a_range = linspace(1e-6,1e-2,1000);  % bubble radius (m)

c = 1485; % speed of sound (m/s)
rho_liq = 1025; % density of liquid (kg/m^3) [water]
rho_air = 1.21;
P_atm = 101.325e3; % atmospheric pressure
g = 9.81; % gravitational acceleration (m/s^2)
d = 220; % water depth (m)

Mm = 16.04e-3; % molar mass of the gas (methane) (kg/mol)
tau0 = 74.5e-3; % surface tension of the gas bubbles (N/m)
tau = 74.5e-3;

mu_liq = 1.519e-3; %shear viscosity N/(m*s)
gamma = 1.299; % heat ratio

T = 281.29; % temperature (K) 
p_v = 872; % vapor pressure of water (Pa)
R = 8.31; %gas constant (m^2.kg.s^-2.K^-1.mol^-1)
C_p = 2.191; % specific heat capacity at const pressure (kJ/kg.K)
K_gas = 8e-2; % thermal conductivity of the gas (W/mK) 
TS = zeros(length(f_range),length(a_range));

for ff = 1:length(f_range)
for aa = 1:length(a_range)

f = f_range(ff);
a = a_range(aa);
w = 2*pi*f; % angular frequency (rad/s)
k = w/c; % wavenumber (1/m)
% k*a
% Thuraisingham model k*a << 1

P_gas = P_atm + rho_liq * g * d + 2*tau/a - p_v; % gas pressure in bubble (Pa) 
rho_gas = Mm / (R*T) * P_gas; % gas density
D_p = K_gas/(rho_gas*C_p); % thermal diffusivity

X = a*sqrt(2*w/D_p); % 
G = gamma / (1 - ((1+1i)*X/2 /(tanh((1+1i)*X/2)) -1)) / (6*1i*(gamma - 1)/X^2);
beta_th =  2*mu_liq / (rho_liq*a^2);% thermal damping (N*s/m)
beta_vis = 3*P_gas / (2*rho_liq*a^2*w) * imag(G); % viscous damping (N*s/m)
beta0 = beta_th + beta_vis;
w0=3 * real(G) * P_gas / (rho_liq*a^2) - 2 * tau0/(rho_liq*a^3);
bs_freq = ((w0/w)^2 - 1 - 2*k*a*beta0/w)^2 + (2*beta0/w + k*a*(w0/w)^2)^2;

% Backscattering cross-section
sigma_bs = (a^2/(bs_freq)) * ((sin(k * a)/k*a)^2/(1+(k*a)^2));

% Target strength
TS(ff, aa) = 10*log10(sigma_bs); %dB re 1 m^2
end
end
%% Plot imagesc: freq x bubble radius x TS
figure;
imagesc(f_range,a_range,TS);
set(gca, 'YDir', 'normal');
colormap("jet")
colorbar;
% title('Range a=[0.01e-4, 0.02e-2]')
% xlim([0.01e-4 0.02e-3])
% saveas(gca, "imagesc_Freq-Radius-TS-1",'png')
% clim([-400 -40])
%% Plot ka x TS
figure;
ka = f_range'*a_range;
kk = 500; % at specific frequency
figure(2)
semilogx(ka(kk,:), TS(kk,:));
%xlim([1e-2 1e2])
xlabel('ka');ylabel('TS (dB re 1 m^2)')
% saveas(gca, "plot_ka-TS",'png')
