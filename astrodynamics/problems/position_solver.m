function F = position_solver(x)

% X = [v1x v1y v1z v2x v2y v2z]

r1 = [1.353683421048372E+08 6.236351063675036E+07  1.616180081151426E+03];
r2 = [-9.694952181127319E+07 4.610774545357175E+07 6.225620771722047E+06];

global f
global g
global f_prima
global g_prima
global p

F(1) = r2(1) - f*r1(1)+g*x(1);
F(2) = r2(2) - f*r1(2)+g*x(2);
F(3) = r2(3) - f*r1(3)+g*x(3);
F(4) = x(4)-f_prima*r1(1)-g_prima*x(1);
F(5) = x(5)-f_prima*r1(2)-g_prima*x(2);
F(6) = x(6)-f_prima*r1(3)-g_prima*x(3);

end