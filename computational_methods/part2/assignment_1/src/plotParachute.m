function plotParachute(time,u,eps,sig,x,Tnod,Tmat)
% Inputs:
%  time  time vector [1 x Ntime]
%  u     nodal displacements matrix [Ndofs x Ntime]
%        (each column corresponds at each timestep)
%  eps   strains matrix [Nelement x Ntime]
%        (each column corresponds at each timestep)
%  sig   stress matrix [Nelemen x Ntime]
%        (each column corresponds at each timestep)
%  x     nodal coordinates matrix [Nnodes x Ndim]
%  Tnod  nodal connectivities matrix [Nelements x NnodesXelement]
%  Tmat  material connectivities matrix [Nelements x 1]
%        It is assumed that mat=1 is for cables and mat=2 for bars


% Get dimensions
Ntime = length(time);
Ndim = size(x,2);
Nnodes = size(x,1);
u = permute(reshape(u,Ndim,Nnodes,Ntime),[2,1,3]);
href = u(1,3,:);
u(:,3,:) = bsxfun(@minus,u(:,3,:),href);

% Conectivities for bars and cables
Tbar   = Tnod(Tmat==2,:);
Tcable = Tnod(Tmat==1,:);

% Global parameters
scale = 1;
current_t = 1;
playBool = false;

for i = 1:Ndim
    Xbar{i} = reshape(x(Tbar,i),size(Tbar))';
    Xcable{i} = reshape(x(Tcable,i),size(Tcable))';
end
Dbar = reshape(sqrt(sum(u(Tbar,:,:).^2,2)),size(Tbar,1),size(Tbar,2),Ntime);
Dcable = reshape(sqrt(sum(u(Tcable,:,:).^2,2)),size(Tcable,1),size(Tcable,2),Ntime);
Ebar = reshape(repelem(eps(Tmat==2,:),size(Tbar,2),1),size(Tbar,2),size(Tbar,1),Ntime);
Ecable = reshape(repelem(eps(Tmat==1,:),size(Tcable,2),1),size(Tcable,2),size(Tcable,1),Ntime);
Sbar = reshape(repelem(sig(Tmat==2,:),size(Tbar,2),1),size(Tbar,2),size(Tbar,1),Ntime);
Scable = reshape(repelem(sig(Tmat==1,:),size(Tcable,2),1),size(Tcable,2),size(Tcable,1),Ntime);

figure('visible','off','color','w','Name','Parachute','position',[50,50,565,540],'SizeChangedFcn',@resizeWindow);

% Create uicontrols
back_btn = uicontrol('style','pushbutton','string','|<','position',[20,20,30,30],'callback',@stepBack);
frwd_btn = uicontrol('style','pushbutton','string','>|','position',[110,20,30,30],'callback',@stepForward);
play_btn = uicontrol('style','pushbutton','string','> / ||','position',[55,20,50,30],'callback',@playPause);
time_sld = uicontrol('style','slider','min',1,'max',Ntime,'Value',1,'position',[145,25,400,20],'callback',@slideTime);
time_lab = uicontrol('style','text','string','time = 0.000 s','position',[20,495,100,20],'backgroundcolor','w');
height_lab = uicontrol('style','text','string','drop = 0.000 m','position',[120,495,100,20],'backgroundcolor','w');
scale_lab = uicontrol('style','text','string','scale = ','position',[240,495,40,20],'backgroundcolor','w');
scale_edt = uicontrol('style','edit','string','1','position',[280,498,40,22],'callback',@scaleFactor);
show_lab = uicontrol('style','text','string','Show:','position',[340,495,40,20],'backgroundcolor','w');
show_pop = uicontrol('style','popupmenu','string',{'None','Displacements','Strain','Stress'},'position',[380,500,160,20],'callback',@showResults);
para_ax = axes('units','pixels','position',[20,70,400,400],'xcolor','none','ycolor','none','zcolor','none');

% Initial plot
p_bar = patch(Xbar{:},zeros(size(Xbar{1})),'edgecolor',[0.5,0.5,0.5],'linewidth',2);
p_cable = patch(Xcable{:},zeros(size(Xcable{1})),'edgecolor',[0.0,0.0,0.0],'linewidth',0.5);
view(30,25);
axis equal;
colormap jet; 
cbar = colorbar;
set(cbar,'visible','off');
set(gcf,'visible','on');

function showResults(varargin)
    updateResults(current_t);
end

function resizeWindow(varargin)
    pos = get(gcf,'position');
    set(time_sld,'position',[145,25,pos(3)-165,20]);
    set(para_ax,'position',[20,70,pos(3)-165,pos(4)-140]);
    set(time_lab,'position',[20,pos(4)-45,100,20]);
    set(height_lab,'position',[120,pos(4)-45,100,20]);
    set(scale_lab,'position',[240,pos(4)-45,40,20]);
    set(scale_edt,'position',[280,pos(4)-42,40,22]);
    set(show_lab,'position',[340,pos(4)-45,40,20]);
    set(show_pop,'position',[380,pos(4)-40,160,20]);
end

function updateTime(tstep)
    set(p_bar,'XData',Xbar{1}+scale*reshape(u(Tbar,1,tstep),size(Tbar))',...
              'YData',Xbar{2}+scale*reshape(u(Tbar,2,tstep),size(Tbar))',...
              'ZData',Xbar{3}+scale*reshape(u(Tbar,3,tstep),size(Tbar))');
    set(p_cable,'XData',Xcable{1}+scale*reshape(u(Tcable,1,tstep),size(Tcable))',...
                'YData',Xcable{2}+scale*reshape(u(Tcable,2,tstep),size(Tcable))',...
                'ZData',Xcable{3}+scale*reshape(u(Tcable,3,tstep),size(Tcable))');
    updateResults(tstep);
    set(time_lab,'string',sprintf('time = %.3f s',time(tstep)));
    set(height_lab,'string',sprintf('drop = %.3f m',href(tstep)));
    drawnow;
end

function updateResults(tstep)
    switch get(show_pop,'value')
        case 1
            set(p_bar,'Cdata',zeros(size(Xbar{1})),'edgecolor',[0.5,0.5,0.5],'linewidth',2);
            set(p_cable,'CData',zeros(size(Xcable{1})),'edgecolor',[0.0,0.0,0.0],'linewidth',0.5);
            set(cbar,'visible','off');
        case 2
            set(p_bar,'Cdata',Dbar(:,:,tstep)','edgecolor','interp','linewidth',2);
            set(p_cable,'CData',Dcable(:,:,tstep)','edgecolor','interp','linewidth',0.5);
            caxis([min([Dbar(:);Dcable(:)]),max([Dbar(:);Dcable(:)])]); 
            set(cbar,'Ticks',linspace(min([Dbar(:);Dcable(:)]),max([Dbar(:);Dcable(:)]),5))
            set(cbar,'visible','on');
        case 3
            set(p_bar,'Cdata',Ebar(:,:,tstep),'edgecolor','flat','linewidth',2);
            set(p_cable,'CData',Ecable(:,:,tstep),'edgecolor','flat','linewidth',0.5);
            caxis([min([Ebar(:);Ecable(:)]),max([Ebar(:);Ecable(:)])]); 
            set(cbar,'Ticks',linspace(min([Ebar(:);Ecable(:)]),max([Ebar(:);Ecable(:)]),5))
            set(cbar,'visible','on');
        case 4
            set(p_bar,'Cdata',Sbar(:,:,tstep),'edgecolor','flat','linewidth',2);
            set(p_cable,'CData',Scable(:,:,tstep),'edgecolor','flat','linewidth',0.5);
            caxis([min([Sbar(:);Scable(:)]),max([Sbar(:);Scable(:)])]); 
            set(cbar,'Ticks',linspace(min([Sbar(:);Scable(:)]),max([Sbar(:);Scable(:)]),5))
            set(cbar,'visible','on');
    end
end

function scaleFactor(varargin)
    scale = str2double(get(scale_edt,'String'));
end

function stepBack(varargin)
    current_t = 1;
    set(time_sld,'value',current_t);
    updateTime(get(time_sld,'value'))
end

function stepForward(varargin)
    current_t = Ntime;
    set(time_sld,'value',current_t);
    updateTime(get(time_sld,'value'))
end

function playPause(varargin)
    if playBool
        playBool = false;
    else
        playBool = true;
        while current_t <= Ntime && playBool
            set(time_sld,'value',current_t);
            updateTime(get(time_sld,'value'))
            current_t = current_t + 1;
            playBool = true;
            pause(time(end)/Ntime);
        end
        playBool = false;
    end
end

function slideTime(varargin)
    current_t = round(get(time_sld,'value'));
    updateTime(current_t);
end

end