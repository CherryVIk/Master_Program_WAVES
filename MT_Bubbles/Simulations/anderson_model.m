%%% Anderson model implementation
%%% Date: 11.03.2024
%%% Author: Viktoriia Boichenko

% Computes an analytical solution for the scattered pressure 
% from an incident plane wave impinging upon a fluid sphere following
% Anderson, V. C., "Sound scattering from a fluid sphere", 
% J. Acoust. Soc. America, 22 (4), pp 426-431, July 1950

f_range = linspace(10e3,1e6,1000); % echosounder freq (Hz=1/s)
a_range = linspace(1e-6,30e-3,1000);  % Oscillations, bubble radius (m)

c_w = 1500; % speed of sound in water (m/s)
c_b = 1.0025*1500; % speed of sound inside bubble (m/s)
rho_liq = 1025; % density of liquid (kg/m^3) [water]
rho_b = 0.66 * rho_liq;
rho_c = rho_liq * c_w;
Theta = 1.571;


% g = 9.81; % gravitational acceleration (m/s^2)
% d = 220; % water depth (m)


TS = zeros(length(f_range),length(a_range));
sigma_bs = zeros(length(f_range),length(a_range));

% for ff = 1:length(f_range)
% for aa = 1:length(a_range)

f = 200 * 1000; % f_range(ff); % acoustic frequency (Hz = s^-1)
a = 3e-3; %a_range(aa); % bubble radius (m)
w = 2*pi*f; % angular frequency (rad/s)
K_ext = w/c_w; % wavenumber in water (1/m)
K_int = w/c_b; % wavenumber in the bubble (1/m)
lambda_w = c_w / f;

Range = 10;
kr = K_ext * Range;
ka = K_ext * a;
ka_p = K_int * a;

g = rho_b / rho_liq;
h = c_b / c_w;
% g = 1.0025; % density contrast
% h = 1.0025; % sound velocity contrast

%% Calculation
m=0;
cosTheta = cos(Theta);

% compute spherical bessel fucntions with ka, k'a, kr, k'r
% Bessel functions 1st degree
jm_kr = sphbes(m, kr);
jm_ka = sphbes(m, ka);
jm_ka_p = sphbes(m, ka_p);

% Bessel functions 2nd degree
nm_kr = sphbes2(m, kr);
nm_ka = sphbes2(m, ka);

% 1st degree derivatives
alpham_ka_p = (2*m+1) * dbesselj(m, ka_p);
alpham_kr = (2*m+1) * dbesselj(m, kr);
alpham_ka = (2*m+1) * dbesselj(m, ka);

% 2nd degree derivatives
betam_kr = (2*m+1) * dbesselj2(m, kr);
betam_ka = (2*m+1) * dbesselj2(m, ka);

% Compute Legendre polynom
Pm = Legendre_pol(m, cosTheta);

% Calculate Cm (no Theta)
Cm_1 = (alpham_ka_p / alpham_ka) * (nm_ka / jm_ka_p) - (betam_ka / alpham_ka) * g*h;
Cm_2 = (alpham_ka_p / alpham_ka) * (jm_ka / jm_ka_p) - g*h;
Cm = Cm_1 / Cm_2;

% Calculate Am and Bm
Am = -(-1i)^m * (2*m+1) / (1 + 1i*Cm);

% Calculate m-value scattered pressure distribution
Ps_n = -Am * Pm * ( jm_kr + 1i*nm_kr);
Pi_n = (-1i)^m * (2*m+1) * Pm * jm_kr;
Ts_n = (-1)^m * (2*m+1) / (1 + 1i*Cm);

% exterior domain particle velocity
%         ThisVel = (-1j / rho_c) * (Am / (2 * m + 1)) * Pm * (Alpham_kr +\
%                   1j * Betam_kr)
%         ThisVel_Inc = (-1j / rho_c) * ((-1j)**m) * Pm * Alpham_kr
%         if m == 0:
%             Ps = ThisPs
%             Pi = ThisPi
%             TS = ThisTS
%             Vel = ThisVel
%             Vel_Inc = ThisVel_Inc
%         else:
%             Ps = Ps + ThisPs
%             Pi = Pi + ThisPi
%             TS = TS + ThisTS
%             Vel = Vel + ThisVel
%             Vel_Inc = Vel_Inc + ThisVel_Inc
%         m += 1
%         Count += 1
%         if max(abs(ThisPs)) < ConvLim and m > 10:
%               Done = True;
%         elif Count > MaxCount:
%               Done = True
%               print('Convergence failure in FluidSphereScat')
% 
%         IsBad = np.isnan(Ps)
%         Ps[IsBad] = 0
% 
%     #Calculate final TS after summation is complete
%     TS = (2 / ka)**2 * abs(TS)**2 * (np.pi * Radius * Radius)
%     TS = 10 * np.log10(TS / ( 4 * np.pi))
%% Functions
function js = sphbes(nu, x)
    % returns the spherical Bessel functions jnu(x)
    % x is a vector or it may be a matrix if nu is a scalar
    % if nu is a row and x a column vector, the output js is a matrix
    [nnu lnu] = size(nu);
    xm = repmat(x, 1, lnu);
    js = sqrt(pi ./(2* xm)) .* besselj(nu + 0.5, x);
end
function js = sphbes2(nu, x)
    % returns the spherical Bessel functions jnu(x)
    % x is a vector or it may be a matrix if nu is a scalar
    % if nu is a row and x a column vector, the output js is a matrix
    [~, lnu] = size(nu);
    xm = repmat(x, 1, lnu);
    js = sqrt(pi ./(2* xm)) .* bessely(nu + 0.5, x);
end
function dJndx = dbesselj(n,x)
    dJndx = n*sphbes(n,x)./x - sphbes(n+1,x);
end
function dJndx = dbesselj2(n,x)
    dJndx = n*sphbes2(n,x)./x - sphbes2(n+1,x);
end
function pol = Legendre_pol(coeffs, degree)
    pol = ones(size(coeffs,1));
    if degree > 0
        pol1 = pol * coeffs;
        if degree == 1
            pol = pol1;
        else
            for ii=2:degree
                pol2 = (coeffs * (2*ii - 1) * pol1 - (ii-1) * pol) / ii;
                pol = pol2;
                pol1 = pol2;
            end
        end
    end
end
%% Plot imagesc: freq x bubble radius x TS
% figure;
% imagesc(f_range,a_range,TS);
% set(gca, 'YDir', 'normal');
% colormap("jet")
% colorbar;
%% Plot ka x TS
% figure(56);
% ka = f_range'*a_range;
% kk = 1; % 30, at specific radius/freq
% semilogx(ka(:,kk), TS(:,kk));
% xlabel('ka');ylabel('TS (dB re 1 m^2)')
% title("Inside function:  ka x TS")