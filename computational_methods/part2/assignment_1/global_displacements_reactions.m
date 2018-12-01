function [u,r] = global_displacements_reactions(KG,f,vr)
% Calculates the global displacements and reaction vectors

vl = setdiff(1:42,vr);

KLL = KG(vl,vl);
KLR = KG(vl,vr);
KRL = KG(vr,vl);
KRR = KG(vr,vr);

FL = f(vl);
FR = f(vr);

ur = zeros(length(vr),1);

ul = KLL \ (FL - (KLR * ur));
r = KRR * ur + KRL * ul - FR;

u(vr,1)=ur;
u(vl,1)=ul;

end

