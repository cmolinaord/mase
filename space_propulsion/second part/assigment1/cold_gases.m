%% Space propulsion
%  Assigment 1 (second part)
%  Carlos Molina

clc; clear all;

%% Input data
g = 9.81; % (m/s^2)
Tc = 311.11; % (K)
Ae_At = [50 80 100];

% Gas characteristics
gas = {'Helium', 'Nitrogen', 'Freon-14', 'Freon-12', 'Ammonia', 'Hydrogen', 'Nitrous oxide'};
gamma = [1.659; 1.40; 1.22; 1.14; 1.31; 1.40; 1.27]; % Gamma constant
M = [4e-3; 28e-3; 88e-3; 121e-3; 17e-3; 2e-3; 44e-3]; % Molar mass
R = [2.0771; 0.2968; 0.0947; 0.06887; 0.4885; 4.1267; 0.18777;]; % Gas constant 


%% Calculations

Pe_P0 = zeros(length(gamma), length(Ae_At));
Ve = zeros(size(Pe_P0));
Isp = zeros(size(Pe_P0));

c_star = ceestar(R, Tc, gamma);

for i = 1:length(gamma)
	for j=1:length(Ae_At)
		% Pressure ratio
		Pe_P0(i,j) = pressure_ratio(gamma(i), Ae_At(j));
		% Exit velocity
		Ve(i,j) = exit_velocity(gamma(i), M(i), Tc, Pe_P0(i,j));
		
		% Specific impulse
		Isp(i,j) = (Ve(i,j)/g) + Ae_At(j)*Pe_P0(i,j)*c_star(i)/g;
	end
end

%% Print results
fprintf('Area ratio:\t')
for j=1:length(Ae_At)
	fprintf('%3.1f\t', Ae_At(j))
end
fprintf('\n');
fprintf('======================================\n')
for i = 1:length(gamma)
	fprintf('%s: \t',cell2mat(gas(i)))
	for j=1:length(Ae_At)
		fprintf('%3.2f\t', Isp(i,j));
	end
	fprintf('\n')
end