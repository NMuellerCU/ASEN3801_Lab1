% Contributors: Nathaniel Mueller
% Course number: ASEN 3801
% File name: main.m
% Created: 1/13/2026

%housekeeping
clear; clc; close all;



%givens
m = 0.05; % kg
Cd = 0.6; % coefficient of drag
g = 9.81; % acceleration due to gravity (m/s^2)
diameter = 0.02; % (m) 
A = pi*(diameter/2)^2; % cross-sectional area (m^2)
h = 1655; %altitude (m)
% 
% rho = stdatmo(h); % density (kg/m^3)
h_vec = [h - 1600, h - 800, h, h + 800, h + 1600]' ; %varying altitude (m)
rho = stdatmo(h_vec)' ;

%time
tspan = [0 20]; % range of time were observing (s)


%initial conditions
p_0 = [0 0 0]' ; % (m)
v_0 = [0 20 -20]' ; % (m/s)

%tolerances
T_a = 1e-8;
T_r = 1e-8;

%initial condions and velocities
x0 = [p_0; v_0 ];
% wind_vel = [0 0 0; 5 0 0;10 0 0; 15 0 0;20 0 0]'; % x
% wind_vel = [0 0 0; 0 5 0;0 10 0; 0 15 0;0 20 0]'; % y
% wind_vel = [ 0 0 0; 0 0 5; 0 0 10; 0 0 15; 0 0 20]'; %  z
wind_vel = [0 0 0]' ; % default

%ODE stuff
opts = odeset('events', @hitGroundEvent, 'RelTol', T_r, 'AbsTol', T_a);
figure;
hold on;
% for wind_conditions = 1:size(wind_vel,2)
%     [t, x] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind_vel(:,wind_conditions)), tspan, x0, opts);

for h_idx = 1:size(h_vec,1)
    [t, x] = ode45(@(t,x) objectEOM(t,x,rho(h_idx),Cd,A,m,g,wind_vel), tspan, x0, opts);


    %plotting
    plot3(x(:,1), x(:,2), -x(:,3), 'LineWidth',1.5);
    hold on;
end
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); 
grid("on")
title('Trajectory of Sphere in Air at 1655m Altitude')
% legend('w = 0i m/s', 'w = 5i m/s', 'w = 10i m/s', 'w = 15i m/s', 'w = 20i m/s');
% exportgraphics(figure(1), 'vary_wind_x.png', 'Resolution', 300);

% legend('w = 0j m/s', 'w = 5j m/s', 'w = 10j m/s', 'w = 15j m/s', 'w = 20j m/s');
% exportgraphics(figure(1), 'vary_wind_y.png', 'Resolution', 300);

% legend('w = 0k m/s', 'w = -5k m/s', 'w = -10k m/s', 'w = -15k m/s', 'w = -20k m/s');
% exportgraphics(figure(1), 'vary_wind_z.png', 'Resolution', 300);

% legend('w = 0k m/s', 'w = -5k m/s', 'w = -10k m/s', 'w = -15k m/s', 'w = -20k m/s');
% exportgraphics(figure(1), 'vary_wind_z.png', 'Resolution', 300);
legend('w = 0k m/s', 'w = 10k m/s', 'w = -10k m/s', 'w = -15k m/s', 'w = -20k m/s');
exportgraphics(figure(1), 'vary_altitude.png', 'Resolution', 300);
%hit ground event
function [position, isterminal, direction] = hitGroundEvent(t, x)
  position = x(3); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 1;   % only detect if going downwards (ie its okay if were going up from below height h to above height h)
end