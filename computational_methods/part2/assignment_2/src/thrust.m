function Te = thrust(x)

	d_1 = 1000; % Drag force at front Spar (N/m)

	y1 = x(1,2);
	y3 = x(53,2);

	fun = @(y) d_1 * (1 - (( (y - y1) / (y3 - y1) ).^2));
	D   = integral(fun,y1,y3); % Drag force (N)

	Te  = D; % Total thrust to compensate Drag (N)

end
