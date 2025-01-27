clear all
close all
clc


% DC Motor Parameters
k = 0.09; % Constant of proportionality (Vs/rad or Nm/A)
R = 7.8;  % Resistance (Ohms)
L = 5e-3; % Inductance (Henries)

% Furuta Pendulum Parameters
L1 = 0.278; % Arm length 1 (m)
L2 = 0.300; % Arm length 2 (m)
l1 = 0.150; % Arm mass offset 1 (m)
l2 = 0.148; % Arm mass offset 2 (m)
m1 = 0.300; % Arm additional mass 1 (kg)
m2 = 0.075; % Arm additional mass 2 (kg)
J1 = 2.48e-2; % Arm moment of inertia 1 (kg·m^2)
J2 = 3.86e-3; % Arm moment of inertia 2 (kg·m^2)
b1 = 1e-4; % Arm viscous damping 1 (N·m·s)
b2 = 2.8e-4; % Arm viscous damping 2 (N·m·s)

% Display parameters
disp('DC Motor Parameters:');
disp(['k = ', num2str(k)]);
disp(['R = ', num2str(R)]);
disp(['L = ', num2str(L)]);

disp('Furuta Pendulum Parameters:');
disp(['L1 = ', num2str(L1)]);
disp(['L2 = ', num2str(L2)]);
disp(['l1 = ', num2str(l1)]);
disp(['l2 = ', num2str(l2)]);
disp(['m1 = ', num2str(m1)]);
disp(['m2 = ', num2str(m2)]);
disp(['J1 = ', num2str(J1)]);
disp(['J2 = ', num2str(J2)]);
disp(['b1 = ', num2str(b1)]);
disp(['b2 = ', num2str(b2)]);
