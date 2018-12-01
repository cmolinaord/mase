function [R, l] = element_R_matrix(x, Tnod, e)
% Gives the length of an element e
dx = x(Tnod(e,1),1) - x(Tnod(e,2),1);
dy = x(Tnod(e,1),2) - x(Tnod(e,2),2);
dz = x(Tnod(e,1),3) - x(Tnod(e,2),3);
l = sqrt(dx*dx+dy*dy+dz*dz);
R = 1/l * [	dx dy dz 0  0  0;
		0  0  0  dx dy dz];
