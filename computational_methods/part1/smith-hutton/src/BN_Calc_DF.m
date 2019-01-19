function [bP]=BN_Calc_DF(AN,AS,AW,AE,bP,nx,ny,Phi_prev,beta);

bP(ny,1)=(1-beta)*(AS(ny,1)*Phi_prev(ny-1,1)+AE(ny,1)*Phi_prev(ny,2));  %Upper-left corner
bP(ny,nx)=(1-beta)*(AS(ny,nx)*Phi_prev(ny-1,nx)+AW(ny,nx)*Phi_prev(ny,nx-1));  %Upper-right corner
bP(1,1)=(1-beta)*(AN(1,1)*Phi_prev(2,1)+AE(1,1)*Phi_prev(1,2)); %Lower-left corner
bP(1,nx)=(1-beta)*(AN(1,nx)*Phi_prev(2,nx)+AW(1,nx)*Phi_prev(1,nx-1));  %Lower-right corner

for i=1:nx
    bP(ny,i)=(1-beta)*(AS(1,i)*Phi_prev(ny-1,i)); %Upper outlet
    bP(i,nx)=(1-beta)*(AW(i,nx)*Phi_prev(i,nx-1)); %Right outlet
end

end