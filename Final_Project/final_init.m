%% Final Project - Equivalent SDOF and Optimal Absorber Parameters
clc; clear; close all;

%% Beam Parameters
% Geometric properties
h = 55e-3;  % Beam height (m)
w = 65e-3;  % Beam width (m)
t = 5e-3;   % Flange thickness (m)
L = 1.5;    % Beam length (m)
a = L/2;    % Machine mass location (m)

% Material properties
E = 210e9;  % Young's modulus (Pa)
nu = 0.33;  % Poisson's ratio
rho = 7850; % Density (kg/m^3)

% Mass properties
mm = 10;    % Machine mass (kg)

% Equivalent cross-sectional properties for I-beam
Iy = (w*h^3)/12 - 2*((t*(h/2)^3)/12);  % Moment of inertia
A = (w*h) - 2*(t*(h/2));  % Cross-sectional area

fprintf('Beam parameters initialized.\n');

%% VAC Parameters
M_VAC = 22;  % Actuator constant (Vs/m = N/A)
R_VAC = 5;   % Actuator resistance (Ohms)
L_VAC = 30e-6;  % Actuator inductance (H)

m_VAC = 0.1;  % Absorber mass (kg)
r_VAC = 50;  % Absorber damping coefficient (Ns/m)
k_VAC = 1000; % Absorber stiffness (N/m)

%% Natural Frequency Calculation
% Fixed-Free beam mode shapes and betaL values
betaL = [1.875, 4.694, 7.855, 10.996]; % First few fixed-free betaL values

% Calculate the first three natural frequencies
omega_n = sqrt((E * Iy) / (rho * A * L^4)) * (betaL(1:3)).^2;
fn = omega_n / (2 * pi);  % Convert to Hz
fprintf('First three natural frequencies: %.2f Hz, %.2f Hz, %.2f Hz\n', fn(1), fn(2), fn(3));

%% Mode Shape Calculation
% Mode shape for fixed-free beam:
bl = betaL(1);
alpha_n = (cosh(bl) + cos(bl)) / (sinh(bl) + sin(bl));

x = linspace(0, L, 200);
C_n = 1 / (cosh(bl*a/L) - cos(bl*a/L) - alpha_n * (sinh(bl*a/L) - sin(bl*a/L)));
w_n = C_n * (cosh(bl*x/L) - cos(bl*x/L) - alpha_n * (sinh(bl*x/L) - sin(bl*x/L)));

% Plot mode shape
figure;
plot(x, w_n, 'b', 'LineWidth', 2);
xlabel('Position along beam (m)');
ylabel('Mode Shape');
title('First Mode Shape of Fixed-Free Beam');
grid on;

%% Equivalent SDOF System
syms x_sym
% Symbolic mode shape
w_x_sym = C_n * (cosh(bl*x_sym/L) - cos(bl*x_sym/L) - alpha_n * (sinh(bl*x_sym/L) - sin(bl*x_sym/L)));

% Equivalent mass calculation
T_sym_beam = int(w_x_sym^2, x_sym, [0, L]);  
T_beam = double(T_sym_beam);
m_beam = T_beam * rho * A * L;

% Effective mass (including machine mass)
m_eq = mm + m_beam;
fprintf('Equivalent mass: %.2f kg\n', m_eq);

% Equivalent stiffness calculation
w_xx = diff(diff(w_x_sym));
V_sym_beam = int(w_xx^2, x_sym, [0, L]);
V_beam = double(V_sym_beam);
k_eq = V_beam * E * Iy / (L^3);
fprintf('Equivalent stiffness: %.2f N/m\n', k_eq);

% System natural frequency
omega_eq = sqrt(k_eq / m_eq);
fn_eq = omega_eq / (2 * pi);
fprintf('System natural frequency (with mass): %.2f Hz\n', fn_eq);

%% Den Hartog's Method - Optimal Absorber Design
mt = 0.05 * m_eq;  % Typical absorber mass ratio (adjustable)
r_opt = 1 / sqrt(2);  % Optimal tuning ratio
omega_t = r_opt * omega_eq; % Tuned frequency

% Optimal damping ratio
zeta_opt = (1/2) * (mt / m_eq) / (1 - r_opt^2);
c_opt = 2 * zeta_opt * sqrt(mt * k_eq);

% Optimal spring stiffness
k_t = mt * omega_t^2;

L_eff = 1/sqrt(k_t/mt);

fprintf('Optimal absorber parameters:\n');
fprintf(' - Absorber mass: %.2f kg\n', mt);
fprintf(' - Absorber natural frequency: %.2f Hz\n', omega_t / (2 * pi));
fprintf(' - Optimal damping coefficient: %.2f Ns/m\n', c_opt);
fprintf(' - Optimal absorber spring stiffness: %.2f N/m\n', k_t);

%% PD Controller Design
Omega = sqrt(k_VAC / m_VAC);
s = 1i * Omega;

num = -(L_VAC * s + R_VAC) * (m_VAC * s^2 + r_VAC * s + k_VAC) - M_VAC^2 * s;

K_P = real(num / M_VAC);
K_D = (1 / Omega) * imag(num / M_VAC);

fprintf('PD Controller Gains:\n');
fprintf(' - K_P: %.4f\n', K_P);
fprintf(' - K_D: %.4f\n', K_D);

sim("beam_undamped.slx")
sim("beam_passive.slx")
sim("beam_active.slx")