clearvars
close all
clc

%% Parameters

%A rectangular beam is considered in this example
E = 210 * 1e9; %Pa
b = 50e-3; %in m
h = 5e-3; %in m

%b = 65e-3
%h = 55e-3
%t = 5e-3

Iy = b*h^3 / 12; %in m^4
A = b*h; %in m^2

%Iy = 4.4557e-07;
%A = t*(h-2*t) + 2*b*t

rho = 7850; %in kg/m^3
l = 1; %in m

%Machine parameters
m_m = 1; %in kg

%% Natural frequencies

beta_L = [3.9266, 7.0685, 10.2101, 13.3518]';
fn = sqrt((E*Iy)/(rho*A*l^4)) .* beta_L.^2 ./ (2*pi)

%% Mode shapes
%Ritz approach:
%w(x) = A sin(beta x) + B cos(beta x) + C sinh(beta x) + D cosh(beta x)
%Apply the coundary conditions (pinned-free) & simplify the expression to solve
%w(x) = C_n * (sin(beta x) + alpha_n * sinh(beta x))

% Calculation of Constant Cn within the mode shape equation depending on the 
% beam support case scaling the mode Wn(x) to 1 at the center of the beam, 
mode = 1;
bl = beta_L(mode);
alpha_n = sin(bl)/sinh(bl);

x_cen = l/2;
wn_cen = 1 * (sin(bl/l*x_cen) + alpha_n * sinh(bl/l*x_cen));
C_n = 1 / (wn_cen);

x = linspace(0, l, 200);
w_n = C_n * (sin(bl/l*x) + alpha_n * sinh(bl/l*x));
plot(x, w_n);

%% Equivalent SDOF
syms x_sym

%Mode shape as symbolic function
w_x_sym = C_n * (sin(bl/l*x_sym) + alpha_n * sinh(bl/l*x_sym));

T_sym_beam = int(w_x_sym^2, x_sym, [0, l]);
T_beam = double(T_sym_beam);

m_beam = T_beam*rho*A*l;
m_beam_system = m_m + m_beam;

w_xx = diff(diff(w_x_sym));
V_sym_beam = int(w_xx^2, x_sym, [0, l]);
V_beam = double(V_sym_beam);
k_beam = V_beam * E*Iy/(l^3);

%Differ because of m_m
omega_beam = sqrt(k_beam/m_beam)
omega_beam_system = sqrt(k_beam/m_beam_system)

f_beam = omega_beam / 2 / pi
f_beam_system = omega_beam_system / 2 / pi