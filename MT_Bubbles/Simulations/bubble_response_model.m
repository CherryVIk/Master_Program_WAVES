function sigma_bs = bubble_response_model(f_range,a_range, model)

% model_type = [1, 'Thuraisingham';
%     2, 'Anderson';
%     3, 'Medwin'];

% f_range = linspace(0.1,300,1000)*1000;
% a_range = 3e-3; % bubble radius (m)
Range = 1;
rho_w = 1025; % density of liquid (kg/m^3) [water]
Theta = 1.571;
c_w = 1500; % speed of sound in water (m/s)

rhoG0 = 0.66; % atmospheric methan density
P_atm = 101.325e3; % atmospheric pressure
g = 9.81; % gravitational acceleration (m/s^2)
d = 20; % water depth (m)
tau = 74.5e-3; % surface tension of the gas bubbles (N/m)
gamma = 1.299; % heat ratio

Pst=P_atm+rho_w*g*d; % static pressure (Pa)
rho_b = rhoG0 * (1 + 2*tau / (Pst * a_range))*(1 + 0.1 * a_range);
c_b = sqrt(gamma*Pst/rho_b); % speed of sound inside bubble (m/s)

if model == 1
    sigma_bs = thuraisingham_model(f_range,a_range, rho_w, c_w);
elseif model == 2
    sigma_bs = anderson_model(f_range,a_range, Range, rho_w, rho_b, Theta, c_w, c_b);
else
    sigma_bs = 0;
end

end