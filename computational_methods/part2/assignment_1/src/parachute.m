close all
clear
clc
%% Input values & constants

n = 14; % Number of nodes
nel = 41; % number of bars
nd = 3; % Number of dimensions
ndof = n * nd; % Number of degrees of freedom

g0 = 9.81; % Gravity constanst (m/s2)
rho_air = 1.225; % Air density (kg/m3)

m_pl = 120; % Mass of the payload (Kg)
[x, Tnod, Tmat, Tdof, Tdn] = get_input_data();

% Material data
% =======================
% Upper surface material
 rhos = 1500; % kg/m3
 S = 17.5; % m2
 Ts = 0.001; % m
 Cd = 1.25;
 m_s = rhos*S*Ts; % Kg


% Elements material properties
m = [% Young modulus [Pa] - Yield strength [Pa] -  Section area [m^2] - Density [kg/m3]
	200000e6  300e6     (0.75^2)*pi/1000000             1500;
	70000e6   240e6     ((6.8^2)-(5.3^2))*pi/1000000    2300];

m_bars = zeros(nel, 1);
for e = 1:nel
	[R,len] = element_R_matrix(x, Tnod, e);
	% Element mass = density * section area * length
	m_bars(e,1) = m(Tmat(e,1),4) * m(Tmat(e,1),3) * len;
end
% Total mass
m_tot = sum(m_bars) + m_s + m_pl;

% SOLVER
% =======================

K = global_stiffness(x,Tnod,m,Tmat,Tdof);

% Global_force: Drag + Weight
% Weight
W = m_tot*g0;

% Terminal velocity: Drag = Weight
Vt = sqrt(2*W/(rho_air*S*Cd));

% Drag at terminal
D = 0.5*rho_air*(Vt^2)*Cd*S;

 %Force vector
f = zeros(ndof,1);
for i = 1:n % For each node
	%Weight
	mbars1 = sum(m_bars(Tnod(:,1) == i))/2;
	mbars2 = sum(m_bars(Tnod(:,2) == i))/2;

	f(i*3,1) = (mbars1 + mbars2) * g0; % (N)
	if i > 5
		% Surface weight addded and Drag substracted
		f(i*3,1) = f(i*3,1) + (m_s/9)*g0 - (D/9);
	end
end
f(3,1) = f(3,1) + m_pl*g0; % Payload weight added

% Imposed (r) and free (l) degrees of freedom
vr = [1 2 3 28 29]; 	 % fixed DoF

[u,r] = global_displacements_reactions(K,f,vr);
% u is the free displacement vector
% r is the reactions vector

%% Visualization
% ======================
[strain,stress] = strain_stress(x,Tnod,m,Tmat,Tdof,u);
plotParachute([1],u,strain,stress,x,Tnod,Tmat)
