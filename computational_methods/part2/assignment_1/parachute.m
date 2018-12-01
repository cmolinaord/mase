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
F0 = m_pl*g0; % Force
[x, Tnod, Tmat, Tdof, Tdn] = get_input_data();

%% External forces
F = zeros(1,ndof);
F(3) = -F0; % Force of the payload
F([18,21,24,27,30,33,36,39,42]) = F0/9; % Upper aerodynamical forces

% Material data
% =======================
% Upper surface material
 rhos = 1500; % kg/m3
 S = 17.5; % m2
 Ts = 0.001; % m
 Cd = 1.25;
 m_s = rhos*S*Ts; % Kg


% Elements material properties
m = [% Young modulus [MPa] - Yield strength [MPa] -  Section area [m^2] - Density [kg/m3]
	200000  300     (0.75^2)*pi/1000000             1500;
	70000   240     ((6.8^2)-(5.3^2))*pi/1000000    2300];

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
for e = 1:n
	%Weight
	mbars1 = sum(m_bars(Tnod(:,1) == e))/2;
	mbars2 = sum(m_bars(Tnod(:,2) == e))/2;

	f(e*3,1) = (mbars1 + mbars2) * g0; % (N)
	if e > 5
		% Surface weight addded and Drag substracted
		f(e*3,1) = f(e*3,1) + (m_s/9)*g0 - (D/9);
	end
end
f(3,1) = f(3,1) + m_pl*g0; % Payload weight added

% Imposed (r) and free (l) degrees of freedom
vr = 1:3; % Fixed DoF
vl = 4:42; % Free DoF
vr = [28 29]; % fixed points
vl = setdiff(1:ndof,vr); % the other points


[u,r] = global_displacements_reactions(K,f,vr);
% u is the free displacement vector
% r is the reactions vector

%% Visualization
% ======================
[strain,stress] = strain_stress(x,Tnod,m,Tmat,Tdof,u);
plotParachute([1],u,strain,stress,x,Tnod,Tmat)
