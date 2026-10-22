% This script file was developed by Toyin Omojola in October 2026
% based on his thesis work at the University of Bath that was submitted in June 2019.
% The code is for steady-state fixed bed kinetic modelling
% Methanol-to-hydrocarbon conversion is used as an example but it can be
% adapted to examples in heterogeneous catalysis


%% housekeeping
tic
clear; close all; clc

%% globalise parameters - equivalent to function file
global k1fValue

%% define initial and final contact times 
tau_initial = 0; tau_final = 20;

%% estimate parameters
T = 370+273;
% specific - equilibirum constant for methanol equilibration reaction -
% equation obtained from literature 
K_p=exp(-9.76+(3200/T)+(1.07*log(T))-(0.66e-3*T)+(0.49e-7*(T^2))+(6500/(T^2)));

k1fValues = 0.0185;
k1bValues = k1fValues/K_p;
k2Values = 0;
k3Values = 0;

%% make globalised parameter equal to estimated parameter
k1fValue = k1fValues; 

%% initial values of mass fraction at zero contact time 
mfo = [0.25;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.75];

%% solve the differential equation
[Tau, mf] = ode15s('MTHmodel', [tau_initial, tau_final], mfo);

% students can check other solvers such as ode23 or ode45 and compare to
% ode15s

%% read experimental data 
% experimental data should be stored in the same folder as script and
% function files
data=xlsread('Experimentaldata.xlsx');

%% plot figures 

toc
