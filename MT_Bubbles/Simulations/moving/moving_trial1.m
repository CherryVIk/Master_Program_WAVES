% [X,Y] = meshgrid(linspace(-5, 5, 50));
% % bubble starts from 
% fcn = @(x,y,k) (x-1)^2 + (y-1)^2*k;
% v = [1:-0.05:-1;  -1:0.05:1];
% for k1 = 1:2
%     for k2 = v(k1,:)
%         surfc(X, Y, fcn(X,Y,k2))
%         axis([-5  5    -5  5    -30  50])
%         drawnow
%         pause(0.1)
%     end
% end

%% Moving red dot in 2D
 % code
x = 150;
y0 = -10;
v = 2;                   %Velocity of 2
figure(2);
for time=1:100
    y= y0 + v*time;        
    plot(x, y,'r*');         %Plot the red dot
    axis([-1 250 -10 10])
    drawnow();
    pause(0.1)
end

%% Moving circle in 2D
 % code
x = 150;
y0 = -10;
v = 2;                   %Velocity of 2
figure(2);
for time=1:100
    y= y0 + v*time;        
    plot(x, y,'r*');         %Plot the red dot
    axis([-1 250 -10 10])
    drawnow();
    pause(0.1)
end

%% Moving sphere spiral in 3D
[X,Y,Z] = sphere;
x = linspace(0, 5, 100);
z = x.^2;
figure(1)
for k2 = 1
    for k1 = 1:length(x)
        
        surf(X+cos(2*pi*x(k1)), Y+sin(2*pi*x(k1)), Z+z(k1))
        hold on
        plot3(cos(2*pi*x), sin(2*pi*x),z,'r-');
        hold off
        axis([-5  5    -5  5    -5  25])
        axis square
        view([-10  20])
        refreshdata
        drawnow
    end
end
%% Moving sphere up in 3D
[X,Y,Z] = sphere;
z = linspace(0, 5, 100);
x = zeros(1,100);
y = zeros(1,100);
% z = x.^2;
figure(1)
for k2 = 1
    for k1 = 1:length(x)
        
        surf(X+x(k1), Y + y(k1), Z+z(k1))
        hold on
        plot3(x, y,z,'r-');
        hold off
        axis([-5  5    -5  5    -5  5])
        axis square
        view([-10  20])
        refreshdata
        drawnow
    end
end