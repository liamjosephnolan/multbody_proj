clearvars
close all
clc

%% Parameter
g = 9.80665; %in m/s^2
m = 2; %in kg
l = 1; %in m

b = 0.75; %in N*m*s/rad
%b = 0; %in N*m*s/rad

phi_init = deg2rad(10); %in rad


%% Run Simscape file

sim("simple_pendulum.slx")

%%

plot(out.phi.time, out.phi.signals.values, '--')
xlabel("t in s")
ylabel("phi in rad")