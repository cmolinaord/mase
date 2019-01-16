function F = gauss_solver(x)

% X = [y s w delta_E p]

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


F(1) = x(1)-(sqrt(mu*x(5))*delta_t)/(r1_n*r2_n*sin(delta_nu));
F(2) = x(2)-((r1_n+r2_n)/(4*sqrt(r1_n*r2_n)*cos(delta_nu/2))-(1/2));
F(3) = x(3)-(mu*(delta_t^2))/(((2*sqrt(r1_n*r2_n)*cos(delta_nu/2)))^3);
F(4) = (x(1)^2)-x(3)/(x(2)+0.5*(1-cos(x(4)/2)));
F(5) = x(1)-(1+(x(4)*sin(x(4)))/(sin(x(4)/2)^3)*(x(2)+(1-cos(x(4)/2))/2));


end