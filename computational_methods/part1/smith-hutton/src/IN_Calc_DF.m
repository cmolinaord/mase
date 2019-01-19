function [AP, bP]=IN_Calc_DF(AP,AN,AS,AW,AE,bP,dx,dy,timestep,n_x,n_y,rho,Phi_prev,beta);
    for i=2:n_x-1
        for j=2:n_y-1
            AP(j,i)=(rho*(dx*dy)/timestep)+beta*(AN(j,i)+AS(j,i)+AE(j,i)+AW(j,i));
            bP(j,i)=(rho*(dx*dy)*Phi_prev(j,i)/timestep)+(1-beta)*((AW(j,i)*(Phi_prev(j,i-1)-Phi_prev(j,i)))-(AE(j,i)*(Phi_prev(j,i)-Phi_prev(j,i+1)))+(AS(j,i)*(Phi_prev(j-1,i)-Phi_prev(j,i)))-(AN(j,i)*(Phi_prev(j,i)-Phi_prev(j+1,i))));
        end
    end
end