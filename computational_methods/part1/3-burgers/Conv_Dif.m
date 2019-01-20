function [R] = Conv_Dif(u, N, Re, LES, Ck, K)
D = zeros(N,1);
C = zeros(N,1);
for k = 1:N
	D(k,1) = Diffusive(k, u, N, Re, LES, Ck);
	for p = -N + k:N
		q = k - p;
		if q ~= 0 && p ~= 0
			if q < 0
				uq = conj(u(-q));
			else
				uq = u(q);
			end
			if p < 0
				up = conj(u(-p));
			else
				up = u(p);
			end
			C(k) = C(k) + up*q*1i*uq;
		end
	end
end
R(K) = D(K)+C(K);
