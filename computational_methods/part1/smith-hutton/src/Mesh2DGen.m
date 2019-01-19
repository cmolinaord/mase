function [posX, posY, dx, dy, n_x, n_y]=Mesh2DGen(n, problem);

switch problem
    case 'SH'
        n_x=n; %Horizontal axis number of nodes
        n_y=n_x/2; %vertical axis number of nodes
        x_ini=-1; %Initial X value for mesh
        x_final=1; %Final X value for mesh
        y_ini=0; %Initial Y value for mesh
        y_final=1; %Final Y value for mesh
    case 'DF'
        n_x=n;
        n_y=n_x; 
        x_ini=0;
        x_final=1; 
        y_ini=0; 
        y_final=1; 
    otherwise
        error ('Wrong problem')
end

%Per each axis, two frontier nodes
dx=(x_final-x_ini)/(n_x-2);
dy=(y_final-y_ini)/(n_y-2);

%Nodes position is defined
posX=zeros(n_x,1);
posX(1)=x_ini;
posX(n_x)=x_final;
for i=2:(n_x-1)
    posX(i)=(posX(1)+dx/2)+dx*(i-2);
end

posY=zeros(n_y,1);
posY(1)=y_ini;
posY(n_y)=y_final;
for i=2:(n_y-1)
    posY(i)=(posY(1)+dy/2)+dy*(i-2);
end

end