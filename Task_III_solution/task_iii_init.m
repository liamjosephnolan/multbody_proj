load("matlab.mat")

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

% Load System Matrices
A = linsys1.A;
B = linsys1.B;
C = linsys1.C;
D = linsys1.D;


% Choose Desired Closed-Loop Poles
desired_poles = [-3 -5 -10 -12 -8];  % Move poles further left but not too extreme




% Compute State Feedback Gain K
K = place(A, B, desired_poles);  % More stable than acker()

% Compute Closed-Loop System
sys_cl = ss(A - B*K, B, C, D);
poles = eig(A - B*K);  % Get the closed-loop poles
disp('Closed-loop poles:');
disp(poles);
