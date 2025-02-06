%% Vorlesungsfile
% Abgabe 02.03.2025 23:55 Uhr
% Hauptteil 3 Seiten
% Berichtssprache: englisch
% part 3: 20% of grading

%% clean up
clear all
close all
clc

%% 03_03_Multibody_Simulation_Task

% DC-Motor Parameters
k_prop = 0.09; % Constant of proportionality (Vs/rad or Nm/A)
R = 7.8;       % Resistance (Ohm)
L = 5e-3;      % Inductance (H)

% Furuta Pendulum Parameters
L1 = 0.278;    % Arm length L1 (m)
L2 = 0.300;    % Arm length L2 (m)
l1 = 0.150;    % Arm mass offset l1 (m)
l2 = 0.148;    % Arm mass offset l2 (m)
m1 = 0.300;    % Arm additional mass m1 (kg)
m2 = 0.075;    % Arm additional mass m2 (kg)
J1 = 2.48e-2;  % Arm moment of inertia J1 (kg·m^2)
J2 = 3.86e-3;  % Arm moment of inertia J2 (kg·m^2)
b1 = 1e-4;     % Arm viscous damping b1 (N·m·s·rad^-1)
b2 = 2.8e-4;   % Arm viscous damping b2 (N·m·s·rad^-1)

r_sphere = 0.025;

%% Calculate controller
% State space model of the system via linearizer
data = load("linsys.mat");
A = data.linsys1.A;
B = data.linsys1.B;
C = data.linsys1.C;
D = data.linsys1.D;

% Controler Design
Q = eye(5);
R = 0.1;

K_lqr = lqr(A,B,Q,R);

%% Run Simscape file(s)

sim("project_task3.slx")