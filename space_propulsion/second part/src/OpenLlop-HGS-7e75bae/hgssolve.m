function [ xc,yc ] = hgssolve(fun,x0,solver,options)
%***********************************************************************************************************
%* HGS 1.3
%* By Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
% 
% Solver selector for hgs code. 
% Motivation: In some situations, the performance or even the functionality
% of hgsTp and hgsisentropic depend on the solver they use. This function 
% is called by hgsT and hsgisentropic and allows the external user to
% select different solvers, including the in-house developped hgsfzero.
% For any issues with the code see the documentation manual.
%
% Inputs:
%   - fun: Function handle of the equation to solve.
%   - x0: Starting point to iterate.
%   - solver: String with the solver to be used
%   - options (optional): Options of the solver to be used. For Matlab solvers use the
%     OPTIMSET structure to generate the options structure. For custom
%     hgsfzero generate the following structure:
%
%       options = struct('x2',5000,'fchange',2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',[]);
%
%     Where info can remain empty if nothing is to be set in screen.
%
% Output:
%   - xc: Solution point x.
%   - yc: Solution point y=f(xc)
%
% See also FSOLVE, FZERO, OPTIMSET
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% initial parameters if option is emtpy
if isempty(options) && strcmp(solver,'hgsfzero')
    options = struct('x2',5000,'fchange',2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',[]);
elseif  isempty(options) % assuming fsolve or fzero
    options = optimset('Display','none');
end

% Get solver as function handle
solv = str2func(solver);

% Launch solver
[xc,yc,flag] = solv(fun,x0,options);

% If flag is negative then run error check code
if flag < 0, ErrorControl(fun,xc,x0,10), end

end

function ErrorControl(f,xc,x0,n)

% Get interval points
x1 = xc -x0/4;
x2 = xc + x0/4;

% Print Error message
fprintf('Problems... plotting function to be solved \n');

% Compute function values where the error has happened
xv=linspace(x1,x2,n); yv(1:n,1) = 0;
for ii=1:length(xv)
    yv(ii)=f(xv(ii));
end

% Plot function to be solved
 figure('Name','hgssolver error control',...
             'Color','w',...
             'NumberTitle','off');
xlabel('x values','Fontsize',16);
ylabel('function values','Fontsize',16);
grid on;

plot(xv,yv,'Linewidth',1.5);
grid

% Raise error
error('uhhh ? can''t solve the equation, check the plot to find the problem');

end