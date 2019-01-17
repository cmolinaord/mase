%% Assignment 1.2:
%% Estimate the ammonia decomposition factor of the Astrium 1N engine in nominal conditions.
% Space propulsion 2019 (MASE-ESEIAT-UPC)
% Carlos Diez - Carlos Molina


% This script computes the estimated ammonia decomposition factor
% of the Astrium 1N engine with data obtained from the specifications of
% this engine by computing the thrust using the equations of hydrazine and
% the properties of the gases during the combustion, and comparing it with
% the theoretical Thrust given by the manufacturer.

%% Code

clear all
clc

% Initial guess for decmposition factor
x0 = 0.5;

opts = optimset('Display','off');
x = fsolve(@thrust_diff, x0, opts);

fprintf('The estimated decomposition factor of ammonia is: %1.4f\n', x)
