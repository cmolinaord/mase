function [a, bP] = BoundaryNodesSH(a,bP, posX, w, alpha)

	for j = 1:w.ny % vertical sides cavity
		a.P(j,1) = 1;
		bP(j,1) = 1-tanh(alpha);
		a.P(j,w.nx) = 1;
		bP(j,w.nx) = 1-tanh(alpha);
	end

	for i = 1:w.nx % Top side cavity
		a.P(w.ny,i) = 1;
		bP(w.ny,i) = 1-tanh(alpha);
	end

	for i = 1:(w.nx/2)  % Inlet
		a.P(1,i) = 1;
		a.N(1,i) = 0;
		a.S(1,i) = 0;
		a.W(1,i) = 0;
		a.E(1,i) = 0;
		bP(1,i) = 1+tanh(alpha*(2*posX(i)+1));
	end

	for i = (w.nx/2):w.nx  % Outlet
		a.P(1,i) = 1;
		a.N(1,i) = 1;
		a.S(1,i) = 0;
		a.W(1,i) = 0;
		a.E(1,i) = 0;
		bP(1,i) = 0;
	end

	% Up-left corner. No north nor west.
	a.P(w.ny,1) = 1;
	a.N(w.ny,1) = 0;
	a.S(w.ny,1) = 0.5;
	a.W(w.ny,1) = 0;
	a.E(w.ny,1) = 0.5;

	% bottom-left corner. No south nor west.
	a.P(1,1) = 1;
	a.N(1,1) = 0.5;
	a.S(1,1) = 0;
	a.W(1,1) = 0;
	a.E(1,1) = 0.5;

	% up-right corner. No north nor east.
	a.P(w.ny,w.nx) = 1;
	a.N(w.ny,w.nx) = 0;
	a.S(w.ny,w.nx) = 0.5;
	a.W(w.ny,w.nx) = 0.5;
	a.E(w.ny,w.nx) = 0;

	% bottom-right corner. No south nor east.
	a.P(1,w.nx) = 1;
	a.N(1,w.nx) = 0.5;
	a.S(1,w.nx) = 0;
	a.W(1,w.nx) = 0.5;
	a.E(1,w.nx) = 0;
end
