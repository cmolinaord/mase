function [D] = Diffusive(k, u, N, Re, LES, Ck)
	if LES == 0
		D = k^2*(1/Re)*u(k);
	else
		m = 2;
		nu_ast = 1+34.5*exp(-3.03*N/k);
		nu_inf = 0.31*(5-m)/(m+1)*sqrt(3-m)*Ck^(-1.5);
		Ekn = u(N)*conj(u(N));
		nu = nu_inf*sqrt(Ekn/N)*nu_ast;
		D = k^2*(1/Re+nu)*u(k);
	end
end
