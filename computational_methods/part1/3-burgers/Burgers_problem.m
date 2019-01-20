function [E, tcomp, u, K] = Burgers_problem(N, Re, At, delta, scheme, LES, Ck )
% Initial Condition and resolution of Burgers problem

%% Initial Conditions

u0 = zeros(N,1);
K = zeros(N,1);
for k = 1:N
	K(k,1) = k;
	if k == 1
		u0(k) = 1;
	else
		u0(k) = 1/abs(k);
	end
end
t = 1;
u(:,t) = u0(:);
R0(:,1) = Conv_Dif(u(:,t), N, Re, LES, Ck, K);
t = t + 1;
u(:,t) = u(:,t-1);
%% Solution
[E, tcomp, u] = Burgers_solver(u, R0, t, At, delta, N, Re, LES, Ck, K, scheme);
end
