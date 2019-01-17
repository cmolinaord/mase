function Pe_P0 = pressure_ratio(gamma, Ae_At)

% gamma: cp/cv
% P0: chamber pressure (bar)
% Pe: exit pressure (bar)

 fun1 = @(Pe_P0)sqrt(0.5*(gamma-1)).*(2./(gamma+1)).^((gamma+1)./(2*(gamma-1)))./((Pe_P0).^(1./gamma).*sqrt(1-(Pe_P0).^((gamma-1)./gamma))) - Ae_At;
 options = optimset('Tolx', 1e-10,...
	 'Tolfun', 1e-10,...
	 'Display','off');
 Pe_P0 = fsolve(fun1,0.4,options);
 
end