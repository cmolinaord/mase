function [V_x_e,V_x_w,V_y_n,V_y_s]=SH_vel(Xpos,Ypos,n_x,n_y);

V_x_e=zeros(n_y,n_x); %X axis velocity (east face)
V_x_w=zeros(n_y,n_x); %X axis velocity (west face)
V_y_n=zeros(n_y,n_x); %Y axis velocity (north face)
V_y_s=zeros(n_y,n_x); %Y axis velocity (south face)

V_x=zeros(n_y,n_x); %X axis velocity (node P)
V_y=zeros(n_y,n_x); %Y axis velocity (node P)

for i=2:n_x-1
    for j=2:n_y-1
        
        deltaXe=(Xpos(i+1)-Xpos(i));
        deltaXw=(Xpos(i)-Xpos(i-1));
        deltaYn=(Ypos(j+1)-Ypos(j));
        deltaYs=(Ypos(j)-Ypos(j-1));
        
        V_x_e(j,i)=2*Ypos(j)*(1-(Xpos(i)+deltaXe/2)^2); %X axis velocity (east face)
        V_x_w(j,i)=2*Ypos(j)*(1-(Xpos(i)-deltaXw/2)^2); %X axis velocity (west face)
        V_y_n(j,i)=-2*Xpos(i)*(1-(Ypos(j)+deltaYn/2)^2); %Y axis velocity (north face)
        V_y_s(j,i)=-2*Xpos(i)*(1-(Ypos(j)-deltaYs/2)^2); %Y axis velocity (south face)
        V_x(j,i)=2*Ypos(j)*(1-Xpos(i)^2); %X axis velocity (node)
        V_y(j,i)=-2*Xpos(i)*(1-Ypos(j)^2); %Y axis velocity (node)
    end
end
return