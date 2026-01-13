% Contributors: Nathaniel Mueller
% Course number: ASEN 3801
% File name: main.m
% Created: 1/13/2026





%initial conditions
clear; clc; close all;
!
%givens
m = 0.05; % kg
c_d = 0.6;
[rho,~,~,~,~,~] = stdatmo(1655);

%initial conditions
p_0 = [0, 0, 0];
v_0 = [0, 20, -20]; % m/s

%tolerances
t_a = 1e-8;
t_r = 1e-8;

%hit ground event
function [position, isterminal, direction] = hitGroundFcn(t, y)
  position = y(1); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end