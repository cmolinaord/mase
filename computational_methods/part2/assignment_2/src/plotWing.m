function plotWing(x, Tnod, u, uint, fint)
	% Inputs:
	%  x     nodal coordinates matrix [Nnodes x Ndim]
	%  Tnod  nodal connectivities matrix [Nelements x NnodesXelement]
	%  u     global displacements/rotations vector [Ndofs x 1]
	%  uint  matrix with element displacements and rotations in LOCAL axes [12 x Nelements]
	%        uint(1,e) : ux for node 1 of element e in local reference frame
	%        uint(2,e) : uy for node 1 of element e in local reference frame
	%        uint(3,e) : uz for node 1 of element e in local reference frame
	%        uint(4,e) : angle_x for node 1 of element e in local reference frame
	%        uint(5,e) : angle_y for node 1 of element e in local reference frame
	%        uint(6,e) : angle_z for node 1 of element e in local reference frame
	%        uint(7,e) : ux for node 2 of element e in local reference frame
	%        uint(8,e) : uy for node 2 of element e in local reference frame
	%        uint(9,e) : uz for node 2 of element e in local reference frame
	%        uint(10,e): angle_x for node 2 of element e in local reference frame
	%        uint(11,e): angle_y for node 2 of element e in local reference frame
	%        uint(12,e): angle_z for node 2 of element e in local reference frame
	%  fint  matrix with element internal forces in LOCAL axes [12 x Nelements]
	%        fint(1,e) : F for node 1 of element e in local reference frame
	%        fint(2,e) : Qy for node 1 of element e in local reference frame
	%        fint(3,e) : Qz for node 1 of element e in local reference frame
	%        fint(4,e) : T for node 1 of element e in local reference frame
	%        fint(5,e) : My for node 1 of element e in local reference frame
	%        fint(6,e) : Mz for node 1 of element e in local reference frame
	%        fint(7,e) : F for node 2 of element e in local reference frame
	%        fint(8,e) : Qy for node 2 of element e in local reference frame
	%        fint(9,e) : Qz for node 2 of element e in local reference frame
	%        fint(10,e): T for node 2 of element e in local reference frame
	%        fint(11,e): My for node 2 of element e in local reference frame
	%        fint(12,e): Mz for node 2 of element e in local reference frame

	% Get dimensions
	Ndim = size(x,2);
	Nnodes = size(x,1);
	Nelements = size(Tnod,1);
	NnodesXelement = size(Tnod,2);
	NdofsXnode = 6;

	% X,Y,Z data
	X = reshape(x(Tnod',1),NnodesXelement,Nelements);
	Y = reshape(x(Tnod',2),NnodesXelement,Nelements);
	Z = reshape(x(Tnod',3),NnodesXelement,Nelements);

	% Displacements
	u = reshape(u,NdofsXnode,Nnodes);
	Ux = reshape(u(1,Tnod'),NnodesXelement,Nelements);
	Uy = reshape(u(2,Tnod'),NnodesXelement,Nelements);
	Uz = reshape(u(3,Tnod'),NnodesXelement,Nelements);
	D = sqrt(Ux.^2+Uy.^2+Uz.^2);

	% Axial deformation
	ux = uint([1,7],:);

	% Deflection Y
	uy = uint([2,8],:);

	% Deflection Z
	uz = uint([3,9],:);

	% Rotation X
	theta_x = uint([4,10],:);

	% Rotation Y
	theta_y = uint([5,11],:);

	% Rotation Z
	theta_z = uint([6,12],:);

	% Axial force
	F = fint([1,7],:);

	% Shear Y
	Qy = fint([2,8],:);

	% Shear Z
	Qz = fint([3,9],:);

	% Torsion
	T = fint([4,10],:);

	% Bending Moment Y
	My = fint([5,11],:);

	% Bending Moment Z
	Mz = fint([6,12],:);

	% Initialize figure
	figure('visible','off','color','w','Name','Wing','position',[50,50,565,540]);

	% Create uicontrols
	scale_lab = uicontrol('style','text','string','scale = ','position',[20,495,40,20],'backgroundcolor','w');
	scale_edt = uicontrol('style','edit','string','1','position',[60,498,40,22],'callback',@scaleFactor);
	show_lab = uicontrol('style','text','string','Show:','position',[120,495,40,20],'backgroundcolor','w');
	show_pop = uicontrol('style','popupmenu','string',{'Displacements','Axial deformation','Deflection Y','Deflection Z',...
	'Rotation X','Rotation Y','Rotation Z','Axial force','Shear force Y','Shear force Z',...
	'Torsion','Bending moment Y','Bending moment Z'},'position',[160,500,160,20],'callback',@showResults);
	para_ax = axes('units','pixels','position',[20,20,400,450]);

	% Initial plot
	sz = 2;
	ori = [x(1,1),x(1,2),x(1,3)]-0.1;
	hold on;
	plot3(ori(1)+[0,sz],ori(2)+[0,0],ori(3)+[0,0],'r','linewidth',1.5); text(ori(1)+1.1*sz,ori(2),ori(3),'X','color','r');
	plot3(ori(1)+[0,0],ori(2)+[0,sz],ori(3)+[0,0],'g','linewidth',1.5); text(ori(1),ori(2)+1.1*sz,ori(3),'Y','color','g');
	plot3(ori(1)+[0,0],ori(2)+[0,0],ori(3)+[0,sz],'b','linewidth',1.5); text(ori(1),ori(2),ori(3)+1.1*sz,'Z','color','b');
	patch(X,Y,Z,ones(size(X)),'edgecolor',[0.5,0.5,0.5],'linewidth',1);
	p = patch(X+Ux,Y+Uy,Z+Uz,D,'edgecolor','interp','linewidth',2);
	view(30,25);
	axis equal;
	colormap jet;
	cbar = colorbar;
	updateCase(D);
	set(gcf,'visible','on');
	set(gca,'xcolor','none','ycolor','none','zcolor','none');

	function showResults(varargin)
		updateResults;
	end

	function resizeWindow(varargin)
		pos = get(gcf,'position');
		set(para_ax,'position',[20,20,pos(3)-165,pos(4)-90]);
		set(scale_lab,'position',[20,pos(4)-45,40,20]);
		set(scale_edt,'position',[60,pos(4)-42,40,22]);
		set(show_lab,'position',[120,pos(4)-45,40,20]);
		set(show_pop,'position',[160,pos(4)-40,160,20]);
	end

	function updateCase(var)
		set(p,'Cdata',var);
		caxis([min(var(:)),max(var(:))]);
		set(cbar,'Ticks',linspace(min(var(:)),max(var(:)),5));
		set(cbar,'visible','on');
	end

	function updateResults
		switch get(show_pop,'value')
		case 1
			updateCase(D);
		case 2
			updateCase(ux);
		case 3
			updateCase(uy);
		case 4
			updateCase(uz);
		case 5
			updateCase(theta_x);
		case 6
			updateCase(theta_y);
		case 7
			updateCase(theta_z);
		case 8
			updateCase(F);
		case 9
			updateCase(Qy);
		case 10
			updateCase(Qz);
		case 11
			updateCase(T);
		case 12
			updateCase(My);
		case 13
			updateCase(Mz);
		end
	end

	function scaleFactor(varargin)
		scale = str2double(get(scale_edt,'String'));
		set(p,'XData',X+scale*Ux,'YData',Y+scale*Uy,'ZData',Z+scale*Uz);
	end

end
