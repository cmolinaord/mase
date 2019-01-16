clear
clc

mu = 1.32712442099e+20; % Sun GM (m3/s2)


r1 = 1.0; % (AU)
r2 = 1.524; % (AU)

a_t = 1.3; % (AU)
e = 1-r1/a_t
% eccentricity

nu_2 = acos((a_t*(1-e^2)-r2)/(e*r2));
nu_2_deg = rad2deg(nu_2)

% Eccentric anomaly at Mars
E2 = acos((e+cos(nu_2)/(1+e*cos(nu_2))))

a_t = 1.94477232e+11; % Semimajor axis of transfer orbit in km

% Time of flight from Earth to Mars
tof = (mu/a_t^3)^(-1/2)*(E2 - e*sin(E2))
