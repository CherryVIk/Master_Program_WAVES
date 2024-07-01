syms f_range, a_range

c = 1500; % speed of sound (m/s)
rho_liq = 1026; % density of liquid (kg/m^3) [water]
rho_air = 1.21;
P_atm = 101.325e3; % atmospheric pressure
g = 9.8; % gravitational acceleration (m/s^2)
d = 5; % water depth (m)
Pst=P_atm+rho_liq*g*d; % static pressure (Pa)

Mm = 16.0e-3; % molar mass of the gas (methane) (kg/mol)
tau0 = 74e-3; % surface tension of the gas bubbles (N/m)
tau = 74e-3;

mu_liq = 1.4e-3; %shear viscosity N/(m*s)
gamma = 1.299; % heat ratio

T = 273+10; % temperature (K) 
p_v = 872; % vapor pressure of water (Pa)
R = 8.31; %gas constant (m^2.kg.s^-2.K^-1.mol^-1)
C_p = 2191; % specific heat capacity at const pressure (kJ/kg.K)
K_gas = 30.6e-2; % thermal conductivity of the gas (W/mK) 
TS = zeros(length(f_range),length(a_range));
sigma_bs = zeros(length(f_range),length(a_range));

for ff = 1:length(f_range)
for aa = 1:length(a_range)

f = f_range(ff);
a = a_range(aa);
w = 2*pi*f; % angular frequency (rad/s)
k = w/c; % wavenumber (1/m)
%% Calculation according to Li2020
P_gas = Pst + (2*tau/a) - p_v;

rho_gas = (Mm/(R*T))*P_gas;
D_p = K_gas/(rho_gas*C_p);

X = sqrt((2*w)/D_p)*a;

Gamma_num1 = (1+1i) * X/2;
Gamma_denum1 = tanh(Gamma_num1);
Gamma_num2 = (6i * (gamma-1))/X^2;
Gamma = gamma/(1-((Gamma_num1/Gamma_denum1)-1)*Gamma_num2);

omegaSquared= (3/(rho_liq*a^2))* (Gamma*P_gas - ((2*tau)/(3*a)));

w_0 = sqrt(real(omegaSquared));

b_th = imag(omegaSquared)/(2*w);
b_vis = 2*mu_liq/(rho_liq*a^2);

b_0 = b_th + b_vis;

sigma_denom2 = (2*(b_0/w)+(w_0^2/w^2)*k*a)^2;
sigma_denom1=((w_0^2/w^2)-1-(2*(b_0/w)*k*a))^2;

sigma_num1=(sin(k*a)/(k*a))^2;
sigma_num2= 1+(k*a)^2;

sigma_bs(ff, aa) = (a^2/(sigma_denom1+sigma_denom2))*(sigma_num1 / sigma_num2);
TS(ff, aa) = 10*log10(sigma_bs(ff, aa)); %dB re 1 m^2
end
end

solve()