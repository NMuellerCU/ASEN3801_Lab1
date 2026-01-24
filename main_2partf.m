% Contributors: Nathaniel Mueller
% Course number: ASEN 3801
% File name: main.m
% Created: 1/13/2026

%housekeeping
clear; clc; close all;



%givens
m = 0.05; % kg
m_vec = linspace(0.001,0.1,100);
Cd = 0.6; % coefficient of drag
g = 9.81; % acceleration due to gravity (m/s^2)
diameter = 0.02; % (m) 
A = pi*(diameter/2)^2; % cross-sectional area (m^2)
h = 1655; %altitude (m)
% 
rho = stdatmo(h); % density (kg/m^3)
% h_vec = [0, 1600, 3200, 4800, 6400]' ; %varying altitude (m)
% rho = stdatmo(h_vec)' ;

%time
tspan = [0 20]; % range of time were observing (s)


%initial conditions
p_0 = [0; 0; 0] ; % (m)
v_0 = [0; 20; -20] ; % (m/s)
v0_hat = v_0/norm(v_0);

KE = 0.5*m*(norm(v_0)^2);
v0_mag_vec = sqrt((2*KE)./m_vec);
v0_vec = v0_hat*v0_mag_vec;

%tolerances
T_a = 1e-8;
T_r = 1e-8;

%initial condions and velocities
% x0 = [p_0; v_0 ];
x0 = p_0*ones(1,length(m_vec));
x0(4:6,:) = v0_vec;
% wind_vel = [0 0 0; 5 0 0;10 0 0; 15 0 0;20 0 0]'; % x
% wind_vel = [0 0 0; 0 5 0;0 10 0; 0 15 0;0 20 0]'; % y
% wind_vel = [ 0 0 0; 0 0 5; 0 0 10; 0 0 15; 0 0 20]'; %  z

% initialize wind matrix set
wind_vel_vec = zeros(3,101,3);
for i=1:2
    wind_vel_vec(i,:,i) = linspace(-100,100,length(wind_vel_vec(i,:,i)));
end
wind_vel = [0 0 0]' ; % default

%ODE stuff
opts = odeset('events', @hitGroundEvent, 'RelTol', T_r, 'AbsTol', T_a);
figure(1)
hold on;
% for wind_conditions = 1:size(wind_vel,2)
%     [t, x] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind_vel(:,wind_conditions)), tspan, x0, opts);

% for h_idx = 1:size(h_vec,1)
%     [t, x] = ode45(@(t,x) objectEOM(t,x,rho(h_idx),Cd,A,m,g,wind_vel), tspan, x0, opts);
% 
% 
%     %plotting
%     plot3(x(:,1), x(:,2), -x(:,3), 'LineWidth',1.5);
%     hold on;
% end
max_distance = 0*m_vec;
for v0_idx = 1:length(m_vec)
    [t, x] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m_vec(v0_idx),g,wind_vel), tspan, x0(:,v0_idx), opts);
    max_distance(v0_idx) = sqrt((x(end,1).^2)+(x(end,2).^2)+(x(end,3).^2));


    %plotting
    %plot3(x(:,1), x(:,2), -x(:,3), 'LineWidth',1.5);
end
%xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
plot(m_vec*1000,max_distance)
xlabel('Mass [g]')
ylabel('Max Distance [m]')
grid("on")
%view(3)
title('Trajectory of Sphere in Air at Varying Masses')
subtitle('Constant Initial Kinetic Energy')
% legend('w = 0i m/s', 'w = 5i m/s', 'w = 10i m/s', 'w = 15i m/s', 'w = 20i m/s');
% exportgraphics(figure(1), 'vary_wind_x.png', 'Resolution', 300);

% legend('w = 0j m/s', 'w = 5j m/s', 'w = 10j m/s', 'w = 15j m/s', 'w = 20j m/s');
% exportgraphics(figure(1), 'vary_wind_y.png', 'Resolution', 300);

% legend('w = 0k m/s', 'w = -5k m/s', 'w = -10k m/s', 'w = -15k m/s', 'w = -20k m/s');
% exportgraphics(figure(1), 'vary_wind_z.png', 'Resolution', 300);

% legend('w = 0k m/s', 'w = -5k m/s', 'w = -10k m/s', 'w = -15k m/s', 'w = -20k m/s');
% exportgraphics(figure(1), 'vary_wind_z.png', 'Resolution', 300);
% legend('Altitude = 0 m', 'Altitude = 1600 m', 'Altitude = 3200 m', 'Altitude = 4800 m', 'Altitude = 6400 m');
%legend(legend_entry,'Location','northeast')
hold off
print("const_KE_variable_mass.png",'-dpng','-r300')
%hit ground event

% Varying Wind
%Initializing Variables
max_distance_matrix = zeros(length(wind_vel_vec(1,:,1)),length(m_vec));
wind_vel_matrix = zeros(length(wind_vel_vec(1,:,1)),length(m_vec));
m_matrix = zeros(length(wind_vel_vec(1,:,1)),length(m_vec));
ideal_mass = zeros(length(wind_vel_vec(1,:,1)),1);
plot_labels = ["Wind in North Direction","Wind in East Direction"];
% k represents the wind direction index
for k=1:2
    for i=1:length(wind_vel_vec(1,:,1)) % i represents the wind speed index
        for j=1:length(m_vec) % j represents the mass index
            % Simulates the trajectory of the sphere given wind and mass
            % conditions
            [t, x] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m_vec(j),g,wind_vel_vec(:,i,k)), tspan, x0(:,j), opts);
            %calulates the landing distance
            max_distance_matrix(i,j) = sqrt((x(end,1).^2)+(x(end,2).^2)+(x(end,3).^2));
            %generates matrix for surface plot
            wind_vel_matrix(i,j) = wind_vel_vec(k,i,k);
            m_matrix(i,j) = m_vec(j);
        end
        %finds the ideal mass for each wind speed
        [~,I] = max(max_distance_matrix(i,:));
        ideal_mass(i) = m_vec(I);
    end
    %plots for each wind direction
    figure(k+1)
    hold on
    surf(m_matrix*1000,wind_vel_matrix,zeros(size(max_distance_matrix)),max_distance_matrix,'EdgeAlpha',0,'FaceColor','interp')
    plot3(ideal_mass*1000,wind_vel_vec(k,:,k),ones(size(ideal_mass)),'Color','k')
    title('Max Distance Based on Mass and Wind Speed')
    subtitle(plot_labels(k) + " and Constant Initial Kinetic Energy")
    xlabel('Mass [g]')
    ylabel('Wind Speed [m/s]')
    view(2)
    grid("on")
    ax = gca;
    ax.Layer = 'top';
    ax.LineWidth = 1;
    h = colorbar;
    h.Label.String = 'Max Distance [m]';
    legend('','Ideal Mass for Max Distance')
    hold off
    print("const_KE_variable_mass_wind"+num2str(k)+".png",'-dpng','-r300')
end
function [position, isterminal, direction] = hitGroundEvent(~, x)
  position = x(3); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 1;   % only detect if going downwards (ie its okay if were going up from below height h to above height h)
end