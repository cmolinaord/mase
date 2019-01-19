function [AP AN AS AW AE bP]=InnerNodesDF(AP,AN,AS,AW,AE,bP,dx,dy,Xpos,Ypos,n_x,n_y,rhoGamma,rho,v0);

%This program is conceived as to use UDS schemes
  V_x=v0*cos(45*pi/180);
  V_y=v0*sin(45*pi/180);

gamma=(rho/rhoGamma);

%Calculation of the inner coefficients (in UDS)
for i=2:(n_x-1)
    for j=2:(n_y-1)
        %Distance between consecutive nodes
        deltaXe=abs(Xpos(i+1)-Xpos(i));
        deltaXw=abs(Xpos(i)-Xpos(i-1));
        deltaYn=abs(Ypos(j+1)-Ypos(j));
        deltaYs=abs(Ypos(j)-Ypos(j-1));
        
        %aN coefficient
        Dn=(gamma*dx)/deltaYn;
        Fn=(rho*V_y)*dx;
        %Pen=Fn/Dn;
        Pen=1; %To impose upwind
        if (-Fn)>=0
            aux=-Fn;
        else
            aux=0;
        end
        AN(j,i)=Dn*Pen+aux;
        
        %aS coefficient
        Ds=(gamma*dx)/deltaYs;
        Fs=(rho*V_y)*dx;
        %Pes=Fs/Ds;
        Pes=1; %To impose upwind
        if Fs>=0
            aux=Fn;
        else
            aux=0;
        end
        AS(j,i)=Ds*Pes+aux;
        
        %aW coefficient
        Dw=(gamma*dy)/deltaXw;
        Fw=(rho*V_x)*dy;
        %Pew=Fw/Dw;
        Pew=1; %To impose upwind
        if (Fw)>=0
            aux=Fw;
        else
            aux=0;
        end
        AW(j,i)=Dw*Pew+aux;
        
        %aE coefficient
        De=(gamma*dy)/deltaXe;
        Fe=(rho*V_x)*dy;
        %Pee=Fe/De;
        Pee=1; %To impose upwind
        if (-Fe)>=0
            aux=-Fe;
        else
            aux=0;
        end
        AE(j,i)=De*Pee+aux;
        
    end
end

return