% Problem 4
% ============

clc
clear all

%% Input positions
r1_v = [1.353683421048372E+08 6.236351063675036E+07 1.616180081151426E+03];
r2_v = [-9.694952181127319E+07 4.610774545357175E+07 6.225620771722047E+06];
mu = 1.33e11;
time_expected = 115; %days
p = 1.2e8;


r1 = norm(r1_v);
r2 = norm(r2_v);
delta_theta = acos(dot(r1_v,r2_v)/(r1*r2));
k = r1 * r2 * (1 - cos(delta_theta) );
l = r1 + r2;
m = r1*r2*(1 + cos(delta_theta));
p1 = k/(l+sqrt(2*m));
p2 = k/(l-sqrt(2*m));
error_cons = 1; % (days)
error = error_cons+1;

%% Iterate
while (error > error_cons)
a = m*k*p / ((2*m - l^2)*p^2 + 2*k*l*p - k^2);
f = 1 - (r2/p)*(1 - cos(delta_theta));
g = r1*r2*sin(delta_theta)/sqrt(mu*p);
E = acos(1 - r1 / a*(1 - f) );
t = g + sqrt(a^3/mu)*(E - sin(E));
tdays = t/(3600*24)
error = abs(-tdays + time_expected);
p = p + 10000;
end

f_p = -sin(E)*sqrt(mu/a)/(r1*r2);
g_p = 1 - (1-cos(E))*a/r2;
v1_v = (r2_v-f*r1_v)/g
norm(v1_v)
v2_v = f_p*r1_v + g_p*v1_v
norm(v2_v)
