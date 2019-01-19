function [AP, AN, AS, AW, AE, bP]=BoundaryNodesSH(AP,AN,AS,AW,AE,bP,posX,n_x,n_y,alpha);


for j=1:n_y %vertical sides cavity
    AP(j,1)=1;
    bP(j,1)=1-tanh(alpha);
    AP(j,n_x)=1;
    bP(j,n_x)=1-tanh(alpha);
end

for i=1:n_x %Top side cavity
    AP(n_y,i)=1;
    bP(n_y,i)=1-tanh(alpha);
end


for i=1:(n_x/2)  %Inlet
           AP(1,i)=1;
           AN(1,i)=0;
           AS(1,i)=0;
           AW(1,i)=0;
           AE(1,i)=0;
           bP(1,i)=1+tanh(alpha*(2*posX(i)+1));
end

for i=(n_x/2):n_x  %Outlet
           AP(1,i)=1;
           AN(1,i)=1;
           AS(1,i)=0;
           AW(1,i)=0;
           AE(1,i)=0;
           bP(1,i)=0;
end


%Up-left corner. No north nor west.
AP(n_y,1)=1;
AN(n_y,1)=0;
AS(n_y,1)=0.5;
AW(n_y,1)=0;
AE(n_y,1)=0.5;
       
%bottom-left corner. No south nor west.
AP(1,1)=1;
AN(1,1)=0.5;
AS(1,1)=0;
AW(1,1)=0;
AE(1,1)=0.5;
       
%up-right corner. No north nor east.
AP(n_y,n_x)=1;
AN(n_y,n_x)=0;
AS(n_y,n_x)=0.5;
AW(n_y,n_x)=0.5;
AE(n_y,n_x)=0;
       
%bottom-right corner. No south nor east.
AP(1,n_x)=1;
AN(1,n_x)=0.5;
AS(1,n_x)=0;
AW(1,n_x)=0.5;
AE(1,n_x)=0;

return