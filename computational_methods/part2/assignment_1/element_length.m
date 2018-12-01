function d = element_length(x, Tnod, e)
% Gives the length of an element e
x1 = x(Tnod(e,1),1);
y1 = x(Tnod(e,1),2);
z1 = x(Tnod(e,1),3);
x2 = x(Tnod(e,2),1);
y2 = x(Tnod(e,2),2);
z2 = x(Tnod(e,2),3);
d = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);

