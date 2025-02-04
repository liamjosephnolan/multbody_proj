%% Finall Project Init
clc

% Beam parameters initialization script

% Geometric properties
h = 55e-3;  % Height (m)
w = 65e-3;  % Width (m)
t = 5e-3;   % Thickness (m)
L = 1.5;    % Beam length (m)
a = L/2;    % Load and mass location (m)

% Material properties
E = 210e9;  % Young's modulus (Pa)
nu = 0.33;  % Poisson's ratio
rho = 7850; % Material density (kg/m^3)

% Mass properties
mm = 10;    % Machine mass (kg)

% Display parameters
fprintf('Beam parameters initialized.\n');
