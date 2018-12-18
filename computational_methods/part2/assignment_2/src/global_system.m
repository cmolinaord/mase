function [u, F] = global_system(data, KG, FG)

% Fixed degrees of freedom at nodes 1 and 5
vr = [1 2 3 4 5 6 25 26 27 28 29 30];
% Imposed degrees of freedom vr
uR = zeros(length(vr),1);
% Free degrees of freedom
vl = setdiff(1:data.ndof,vr);

% Partitioned system
KLL = KG(vl,vl);
KLR = KG(vl,vr);
KRL = KG(vr,vl);
KRR = KG(vr,vr);
FGL = FG(vl);
FGR = FG(vr);

% System resolution
uL = inv(KLL) * (FGL - KLR * uR);
RR = KRR*uR + KRL*uL - FGR;

% Obtain generalized displacement vector
u(vl,1) = uL;
u(vr,1) = uR;

% Total forces vector %%
Reactions       = zeros(length(FG),1);
Reactions(vr)   = RR;
F               = FG + Reactions;

end
