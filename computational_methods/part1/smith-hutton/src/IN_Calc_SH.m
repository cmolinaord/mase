function [AP, bP]=IN_Calc_SH(AP,AN,AS,AW,AE,bP,dx,dy,timestep,n_x,n_y,rho,Phi_prev,beta);

source=0; %As per problem specifications

%The inner nodes are updated

for i=2:n_x-1
    for j=2:n_y-1
        
        %aP coefficient update
        AP(j,i)=(rho*(dx*dy)/timestep)+beta*(AN(j,i)+AS(j,i)+AE(j,i)+AW(j,i))-source*(dx+dy);
        
        %bP coefficient update
        bP(j,i)=(rho*(dx*dy)*Phi_prev(j,i)/timestep)+(source*(dx+dy))+(1-beta)*((AW(j,i)*(Phi_prev(j,i-1)-Phi_prev(j,i)))-(AE(j,i)*(Phi_prev(j,i)-Phi_prev(j,i+1)))+(AS(j,i)*(Phi_prev(j-1,i)-Phi_prev(j,i)))-(AN(j,i)*(Phi_prev(j,i)-Phi_prev(j+1,i))));
        
    end
end


end