function sigma_bs = bubble_response_model(f_range,a_range, model)

% model_type = [1, 'Thuraisingham';
%     2, 'Anderson';
%     3, 'Medwin'];

c_w = 1500; % speed of sound in water (m/s)
c_b = 1.0025*1500; % speed of sound inside bubble (m/s)
rho_w = 1026; % density of liquid (kg/m^3) [water]
rho_b = 0.66 * rho_w;
Theta = 1.571;
Range = 1;
if model == 1
    sigma_bs = thuraisingham_model(f_range,a_range, Range, rho_w, rho_b, Theta, c_w, c_b);
elseif model == 2
    sigma_bs = anderson_model(f_range,a_range, Range, rho_w, rho_b, Theta, c_w, c_b);
else
    sigma_bs = 0;
end

end