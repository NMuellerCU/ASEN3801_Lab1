clear; clc; close all;

beta = 60*pi/180;
theta_i = beta;
theta_j = beta;
theta_k = beta;
R1 = [1, 0, 0;
            0, cos(theta_i), sin(theta_i);
            0, -sin(theta_i), cos(theta_i)];
R2 = [cos(theta_j), 0, -sin(theta_j);
            0, 1, 0;
            sin(theta_j),0, cos(theta_j)];
R3 = [cos(theta_k), sin(theta_k), 0;
            -sin(theta_k), cos(theta_k), 0;
            0, 0, 1];
% 1 2 3
rot_matrix_1 = R3 * R2 * R1; %123

[V1, D1] = eig(rot_matrix_1);
disp(rot_matrix_1);
disp(V1)
disp(D1)
% 3 1 3
rot_matrix_2 = R3 * R1 * R3; %313

[V2, D2] = eig(rot_matrix_2);
disp(rot_matrix_2);
disp(V2)
disp(D2)
% 3 2 1
rot_matrix_3 = R1 * R2 * R3; %321

[V3, D3] = eig(rot_matrix_3);
disp(rot_matrix_3);
disp(V3)
disp(D3)