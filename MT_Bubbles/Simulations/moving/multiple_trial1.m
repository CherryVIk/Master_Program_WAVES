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
% https://de.mathworks.com/matlabcentral/answers/460216-generate-n-random-uniformly-distributed-points-in-a-specific-area
%% with rand intervals as rectangles around???
figure(30)
a=rand(2,100)*5;
scatter(a(1,:),a(2,:))
%% Generate N random uniformly distributed points in a square  frame
nreq = 1000;
% an extra 20% in case we reject too many.
rejectfrac = (1.5^2 - 1)/4;
oversample = 0.2;
xy = [];
nxy = 0;
center = [0, 0];
side_x = 1;side_y = 2;
side2_x = side_x*1.5; 
side2_y = side_y*1.5;
while nxy < nreq
  % assume we have a square incribed in a square. So the square has edge length of 1.5*side_in.
  nsample = ceil((1 + oversample)*(nreq - nxy)/rejectfrac);
  % this mext line uses a feature found in R2016b or later.
  % rand num in [-1, 1]*side_in
  xypoints = (rand(nsample,2)*2 - 1)*side_x + center;

  % points in the vicinity of the center square are deleted
  if xypoints(:,1) > side_x/2 + center(1) && xypoints(:,1) > center(1) - side_x/2
      if xypoints(:,2) < side_y/2 + center(2) && xypoints(:,1) > center(2) - side_y/2
      end
  end
  % delete points inside the circle.
  xypoints(sum((xypoints - center).^2,2) <= side_x^2,:) = [];
  
  % were there enough points
  xy = [xy;xypoints];
  nxy = size(xy,1);
  if nxy > nreq
    xy = xy(1:nreq,:);
  end
end
figure(10)
plot(xy(:,1),xy(:,2),'.')
axis equal
grid on
yline(0);
xline(0);