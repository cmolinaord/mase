function [Phi_new]=GS_solver(AP,AN,AS,AW,AE,bP,n_x,n_y,rho,Phi_ini,iter_lim,delta,Phi_new,beta);
%Gauss-Seidel

dif=1000;
iter=0;
aux=0;

while iter<iter_lim && dif>delta
    aux=0;
    dif=0;
    iter=iter+1;
    
    for i=2:n_x-1 %Inner nodes
        for j=2:n_y-1
            
            Phi_new(j,i)=(1/AP(j,i))*(beta*(AE(j,i)*Phi_ini(j,i+1)+AW(j,i)*Phi_ini(j,i-1)+AN(j,i)*Phi_ini(j+1,i)+AS(j,i)*Phi_ini(j-1,i))+bP(j,i));    
         
        end
    end
    
    
    for i=2:n_x-1 %Top and bottom, no corners
        %Top:
        Phi_new(n_y,i)=(1/AP(n_y,i))*(beta*(AE(n_y,i)*Phi_ini(n_y,i+1)+AW(n_y,i)*Phi_ini(n_y,i-1)+AS(n_y,i)*Phi_ini(n_y-1,i))+bP(n_y,i));
        
        %Bottom:
        Phi_new(1,i)=(1/AP(1,i))*(beta*(AE(1,i)*Phi_ini(1,i+1)+AW(1,i)*Phi_ini(1,i-1)+AN(1,i)*Phi_ini(1+1,i))+bP(1,i));    
    end
    
    for j=2:n_y-1 %Left and right, no corners
        %Left:
        Phi_new(j,1)=(1/AP(j,1))*(beta*(AE(j,1)*Phi_ini(j,1+1)+AN(j,1)*Phi_ini(j+1,1)+AS(j,1)*Phi_ini(j-1,1))+bP(j,1));    
        
        %Right:
        Phi_new(j,n_x)=(1/AP(j,n_x))*(beta*(AW(j,n_x)*Phi_ini(j,n_x-1)+AN(j,n_x)*Phi_ini(j+1,n_x)+AS(j,n_x)*Phi_ini(j-1,n_x))+bP(j,n_x));    
    end
    
    %Corners:
    %Up-left
    Phi_new(n_y,1)=(1/AP(n_y,1))*(beta*(AE(n_y,1)*Phi_ini(n_y,1+1)+AS(n_y,1)*Phi_ini(n_y-1,1))+bP(n_y,1));    
    %Up-right:
    Phi_new(n_y,n_x)=(1/AP(n_y,n_x))*(beta*(AW(n_y,n_x)*Phi_ini(n_y,n_x-1)+AS(n_y,n_x)*Phi_ini(n_y-1,n_x))+bP(n_y,n_x));    
    %Bottom-left:
    Phi_new(1,1)=(1/AP(1,1))*(beta*(AE(1,1)*Phi_ini(1,1+1)+AN(1,1)*Phi_ini(1+1,1))+bP(1,1));
    %Bottom-right:
    Phi_new(1,n_x)=(1/AP(1,n_x))*(beta*(AW(1,n_x)*Phi_ini(1,n_x-1)+AN(1,n_x)*Phi_ini(1+1,n_x))+bP(1,n_x));    
         
   
    for i=1:n_x
        for j=1:n_y
            aux=abs(Phi_new(j,i)-Phi_ini(j,i));
            if aux>dif
                dif=aux;
            end
            
        end
    end
    
    if dif>delta
        Phi_ini=Phi_new;
    end
    
end

end