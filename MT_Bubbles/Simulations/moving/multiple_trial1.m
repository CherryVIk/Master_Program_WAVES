%%% Data: 29.02.2024
%%% Thema: Create multiple bubbles and show it in 2D
close all
clear 

bubble_seed = 42;
Nbubbles = 10; % number of bubbles
Ndims = 2; % [x y] or [x y z]
bubble_step = 0.1; % one step of bubble moving
bubble_Dist = 1; % maximum distance

% posTar = [0 40;10 20]; % [x y]
%% Initialise data

a=-5;
b=5;
% Generate values from the uniform distribution on the
%        interval (a, b).
% rng(bubble_seed);
x = a + (b-a).*rand(Nbubbles,1);
y = a + (b-a).*rand(Nbubbles,1);
z = zeros(Nbubbles, 1);

figure(1)
plot3(x,y,z, '-ok')
% axis([a b a b a b])
grid on

%% Move points vertically 
figure(2);
hold on
for jj=0:bubble_step:bubble_Dist
%     x = x + jj;
%     y = y + jj;
    z = z + jj;
    plot3(x,y,z, '-ok')
    grid on
    pause(0.05)
    drawnow
end
hold off