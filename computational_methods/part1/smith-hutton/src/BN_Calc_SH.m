function [bP]=BN_Calc_SH(AN,AS,AW,AE,bP,n_x,n_y,Phi_prev,beta);

%% Update of the bottom outlet
for i=(n_x/2):n_x
    bP(1,i)=(1-beta)*(AN(1,i)*Phi_prev(2,i));
end


%% Update of the corners

%Up-left:
bP(n_y,1)=(1-beta)*(AS(n_y,1)*Phi_prev(n_y-1,1)+AE(n_y,1)*Phi_prev(n_y,2));

%Up-right:
bP(n_y,n_x)=(1-beta)*(AS(n_y,n_x)*Phi_prev(n_y-1,n_x)+AW(n_y,n_x)*Phi_prev(n_y,n_x-1));

%Bottom-left:
bP(1,1)=(1-beta)*(AN(1,1)*Phi_prev(2,1)+AE(1,1)*Phi_prev(1,2));

%Bottom-right:
bP(1,n_x)=(1-beta)*(AN(1,n_x)*Phi_prev(2,n_x)+AW(1,n_x)*Phi_prev(1,n_x-1));

end