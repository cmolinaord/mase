% Problem 5

a = -36000e3; % Semimajor axis (m)
e = 1.1823; % Eccentricity

nu_1 = deg2rad(15);		% (rad)
nu_2 = deg2rad(120);	% (rad)

F1 = 2*atanh(tan(0.5*nu_1)*((e+1)/(e-1))^(-0.5))
F2 = 2*atanh(tan(0.5*nu_2)*((e+1)/(e-1))^(-0.5))

mu = 398600441800000.0

t = 1/sqrt(mu/(-a)^3) * ((e*sinh(F2)-F2) - (e*sinh(F1)-F1))
