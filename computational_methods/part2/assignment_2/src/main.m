%%
% Computational engineering
% MASE (Master's Degree in Space and Aeronautical Engineering)
% UPC- ESEIAAT (Terrassa)
% Carlos Molina (carlosmolina.ord@gmail.com) 17-Dec-2018
% https://github.com/cmolinaord/mase/tree/master/computational_methods/part2/assignment_1/src
%

close all;
clear all;
clc;

%% Input data

[x, Tnod, Tmat, dat, Tn, mat, data] = input_wing;

%% Step 1:
% Associated parameters calculation

[L, A, I_y, I_z, J, V, m] = associated_parameters(data, x, Tnod, dat, Tmat, mat);

%% Step 2:
% Mass & density calculations

[rho_hat, rho_eff] = density(data, Tmat, m, V, mat);

%% Step 3:
% Thrust.

Te = thrust(x);
fprintf('The thrust provided by the engine is %2.3fkN\n', Te*1e-3)

%% Step 4:
% Stiffness element and global matrix
[Ke, KG, R, Re] = stiffness_matrix(data, Tmat, Tn, mat, dat, A, L, I_z, I_y, J);

%% Step 5:
% Element and global force vector

[Fe, FG] = force_vector(data, x, Tnod, L, Te, rho_eff, A, R, Re, Tn);

%% Step 6:
% Solve the global system

[u, F] = global_system(data, KG, FG);

%% Step 7:
% Global displacements, rotations and forces for each element

[uint, fint] = local_DisRotFor(data, Tn, F, u, Re, Ke);

%% Step 8:
% Plot final results

plotWing(x,Tnod,u,uint,fint);
