%% Assignment 1.2:
%% Estimate the ammonia decomposition factor of the Astrium 1N engine in nominal conditions.
% Space propulsion 2019 (MASE-ESEIAT-UPC)
% Carlos Diez - Carlos Molina


% This script computes the estimated ammonia decomposition factor
% of the Astrium 1N engine with data obtained from the specifications of
% this engine by computing the thrust using the equations of hydrazine and
% the properties of the gases during the combustion, and comparing it with
% the theoretical Thrust given by the manufacturer.

%% Manufacturer specifications

Isp = 220;			% Nominal specific impulse (s)
m_tot = 0.29;		% Total mass (kg)
deltah_vap = 44.5;	% Vaporization heat of hydrazine (kJ/mol)



%% Code

% Initial guess for decmposition factor
x0 = 0.5;

opts = optimset('Display','off');
x = fsolve(@thrust_diff, x0, opts);

fprintf('The estimated decomposition factor of ammonia is: %1.4f\n', x)