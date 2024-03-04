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
%%  a bivariate normal distribution
mu = [1 2];  %mean vector
Sigma = [1 .5; .5 2]; R = chol(Sigma); %covariance vector
z = repmat(mu,100,1) + randn(100,2)*R;
%% uniform distribution
% Generate integer values from the uniform distribution on the set 
% rng(bubble_seed);
%  constrains a, b
a=30;
b=35;
Nbubbles=10;
posTar = randi([a,b],Nbubbles,3);
% posTar(:,3) = posTar(:,3) + move_ii;
move_ii = 10;
posTar = posTar + move_ii*ones(Nbubbles, 3);
NTargets = size(posTar, 1);
bDirectSound = 0;

x = posTar(:,1);
y = posTar(:,2);
z = posTar(:,3);
figure
plot3(x,y,z, '-ok')
grid on
%% Generate N random uniformly distributed points in a specific area
N=200; %number of points
n=1;   %Iterator
x_range=[-1 1]; %Range of width
mid_point=[mean(x_range),mean(x_range)]; %Center of box
radius=1;   %Radius of circle
point_arr=zeros(N,2); %This will hold the point
while n<=N
    temp_xy = (x_range(2)-x_range(1)).*rand(1,2) + x_range(1); %Generate 2 random numbers x and y
    d = sqrt(sum((temp_xy-mid_point).^2)); %Find distance between the point generated and the middle
    if d>radius %If the distance is smaller than the radius, discard it, else accept it
        point_arr(n,:)=temp_xy;
        n=n+1;
    end
end
scatter (point_arr(:,1),point_arr(:,2))