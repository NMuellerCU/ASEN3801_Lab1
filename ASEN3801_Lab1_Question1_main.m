%% ASEN 3801 Lab 1

% Contributors: Evan Karen
% Course number: ASEN 3801
% File name: ASEN_Lab1_main.m
% Created: 1/13/2026

%Initializing program
 clc;
 close all;
 clear;

 %% Problem 1

 % a.)

 %state vector derivative function call
 dy = @(t,y) [-9*y(1) + y(3);...
     4*y(1)*y(2)*y(3) - (y(2))^2;...
     2*y(1) - y(2) - 2*y(4);...
     y(2)*y(3) - (y(3))^2 - 3*(y(4))^3];

 y0 = [1;1;1;1]; %Initial Condition
 tspan = [0,20]; %Time span for dyanmical system
 names = ["$\bf{w}$","$\bf{x}$","$\bf{y}$","$\bf{z}$"]; %Names for each 
 % variable (names(1) corresponds to y(1))
 opts = odeset('RelTol',10^-(8),'AbsTol',10^-(8)); %Sets Tolerance values 
 % for initial ode45 call

 [t,y] = ode45(@(t,y) dy(t,y),tspan,y0,opts); %Calls ode45 to numerically 
 % solve the dynamical system; returns t, a time vector, and y, a 4 vector 
 % with values corresponding to every t value

 % Generates figure for Problem 1 part a
 figure(1)
 set(gcf, 'Position', [0 50 1000 600]) %Sets the figure's size and position
 for i=1:4 %For each vector component, plots it's values over time on a 
     % different subplot
    subplot(4,1,i)
    hold on
    ylabel(names(i),'Interpreter','latex','Rotation',0,'FontSize',24)
    plot(t,y(:,i),'LineWidth',2,'Color','k')
    grid('on')
 end
 %additional plot formatting
 xlabel('$\bf{t}$','Interpreter','latex','FontSize',24)
 sgtitle('Solution to Dynamical System - Problem 1.a')
 hold off
 print("ASEN3801_Lab1_Problem_1a",'-dpng','-r300')

 % b.)
 opts = odeset('RelTol',10^-(12),'AbsTol',10^-(12)); %Reference Tolerance 
 % values

 [tR,yR] = ode45(@(tR,yR) dy(tR,yR),tspan,y0,opts); %returns the 
 % numerically solved reference dynamical system

 tol = [10^(-2),10^(-4),10^(-6),10^(-8),10^(-10),10^(-12)]; %array of 
 % tolerances

 vals = zeros([4,6]); %initializes solution table as a empty table

 % for loop running through each tolerance value
 for i=1:length(tol)
     opts = odeset('RelTol',tol(i),'AbsTol',tol(i)); %sets tolerance values
     [t,y] = ode45(@(t,y) dy(t,y),tspan,y0,opts); %runs ode45 for 
     % (wR,xR,yR,zR) vector with given tolerance values
     vals(:,i) = abs(y(end,:)' - yR(end,:)'); %finds the difference of the 
     % t=20 values for each vector (w,x,y,z) and tolerance
 end

 