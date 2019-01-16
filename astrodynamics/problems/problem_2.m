
% Julian year. 365.25 days of 24h of 3600s
sec_yr = 31557600;

%% Input data
e = 0.21634; % eccentricity
T = 4.3856; % Period (yr)
T = T*sec_yr; % (s)

dT = 1.2841; % Time delta (yr)
dT = dT*sec_yr; % (s)

% Constants

mu_sun = 1.32712442099e+20; % Solar GM (m3/s2)


%% Computations
% Semimajor-axis
a = ((T/2/pi)^2*mu_sun)^(1/3);

% Mean anomaly equation
fun = @(E) (mu_sun/a^3)^(-0.5)*(E-e*sin(E))-dT;

% Solve to find E eccentric anomaly (result in rad)
options = optimset('Display','iter');
options.TolFun = 1e-8;
options.MaxIter = 200;
E1 = fsolve(fun, 2);

fprintf('%1.5f rad\n%3.5f deg\n', E1, rad2deg(E1))

