% Computational Fluid Dynamics
% Burger's equation
% MASE (Master's Degree in Space and Aeronautical Engineering)
% Carlos Molina
% January 2019

clc
clear
close all

%% Initial data

Re = 40;
N = 80;
alpha = 0.01;
Ck = 0.05;  %Kolmogorov constant
% Ck = 0.4523;

LES = 0;  %0 to deactivate, 1 to activate

timestep = alpha*Re/N^2; % time step
delta_conv = 1e-6;
%% Initial Conditions

u0 = zeros(N,1);
K = zeros(N,1);
for k = 1:N
	K(k,1) = k;
	u0(k) = 1/abs(k);
end
t = 1;
u(:,t) = u0(:);
R_ini(:,1) = Conv_Dif(u(:,t), N, Re, LES, Ck, K);
t = t+1;
u(:,t) = u(:,t-1);

%% Solution
[E, tcomp, u] = Burgers_solver(u, R_ini, t, timestep, delta_conv, N, Re, LES, Ck, K);

%% Post-Process
slope(K) = 1./(K.^2);
figure()
loglog(K, E, '-+', K, slope, '-.r')
xlabel('k')
ylabel('E_k')
grid
if LES == 0
	l_les = '';
else
	l_les = [' LES C_K = ',num2str(Ck)];
end
l_num = ['N = ',num2str(N),l_les];
legend(l_num,'Slope = -2');
