%% LAMBERT'S PROBLEM: GAUSS' METHOD %%

clc
clear all
close all

mu = 1.32712440018e11; %sun, km3/sec
delta_t = 115*24*3600; % time between two dates (October 8 1989 and Feb 10 1990), sec
maxiter = 1000;
error_cons = 1e-15;
error = error_cons+1;


r1 = [1.353683421048372E+08 6.236351063675036E+07  1.616180081151426E+03];
r2 = [-9.694952181127319E+07 4.610774545357175E+07 6.225620771722047E+06];
r1_n = norm(r1);
r2_n = norm(r2);

r1xr2 = cross(r1,r2);
if r1xr2(3) < 0
   delta_nu = 2*pi - acos(r1*r2'/(r1_n*r2_n));
end

if r1xr2(3) >= 0
   delta_nu =  acos (r1*r2'/(r1_n*r2_n));
end

x0 = [100 10 10 0.5 500];
x = fsolve(@gauss_solver, x0);

global f
global g
global f_prima
global g_prima

y = x(1);
s = x(2);
w = x(3);
delta_E = x(4);
p = x(5);

f = 1-(r2_n/p)*(1-cos(delta_nu));
g = (r1_n*r2_n*sin(delta_nu))/sqrt(mu*p);
f_prima = sqrt(mu/p)*tan(delta_nu/2)*(((1-cos(delta_nu))/p)-(1/r1_n)-(1/r2_n));
g_prima = 1-(r1_n/p)*(1-cos(delta_nu));

x02 = [10 10 10 10 10 10];
velocities = fsolve(@position_solver,x02);

v1 = [velocities(1) velocities(2) velocities(3)];
v2 = [velocities(4) velocities(5) velocities(6)];

v1n = norm(v1);
v2n = norm(v2);

