function [ E,tcomp,u ] = Burgers_solver(u, R0, t, At, delta, N, Re, LES, Ck, K)
% Burgers equation solver applying full explicit scheme

tic
conver = 0;
while conver == 0
	R(:,1) = Conv_Dif(u(:,t), N, Re, LES, Ck, K);
	u(K,t+1) = u(K,t)-At*R(K,1);
	u(1,t+1) = 1;
	if max(abs(u(:,t+1) - u(:,t)))<delta
		conver = 1;
		E(K) = u(K,t+1).*conj(u(K,t+1));
	else
		t = t+1;
		R0(K) = R(K);
	end
end
tcomp = toc;
end
