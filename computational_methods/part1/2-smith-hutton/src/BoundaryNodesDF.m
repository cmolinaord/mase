function [a bP] = BoundaryNodesDF(a,bP,w)

for j = 1:w.ny %vertical sides cavity
	a.P(j,1) = 1; %Left
	bP(j,1) = 1;
	a.P(j,w.nx) = 1; %Right
	a.W(j,w.nx) = 1;
	bP(j,w.nx) = 0;
end

for i = 1:w.nx %horizontal side cavity
	a.P(w.ny,i) = 1;%Upper
	a.S(w.ny,i) = 1;
	bP(w.ny,i) = 0;
	a.P(1,i) = 1; %Bottom
	bP(1,i) = 0;
end

%Up-left corner. No north nor west.
a.P(w.ny,1) = 1;
a.N(w.ny,1) = 0;
a.S(w.ny,1) = 0.5;
a.W(w.ny,1) = 0;
a.E(w.ny,1) = 0.5;

%bottom-left corner. No south nor west.
a.P(1,1) = 1;
a.N(1,1) = 0.5;
a.S(1,1) = 0;
a.W(1,1) = 0;
a.E(1,1) = 0.5;

%up-right corner. No north nor east.
a.P(w.ny,w.nx) = 1;
a.N(w.ny,w.nx) = 0;
a.S(w.ny,w.nx) = 0.5;
a.W(w.ny,w.nx) = 0.5;
a.E(w.ny,w.nx) = 0;

%bottom-right corner. No south nor east.
a.P(1,w.nx) = 1;
a.N(1,w.nx) = 0.5;
a.S(1,w.nx) = 0;
a.W(1,w.nx) = 0.5;
a.E(1,w.nx) = 0;
end
