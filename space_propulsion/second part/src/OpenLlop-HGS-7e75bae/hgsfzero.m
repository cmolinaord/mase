function [xc, yc, flag] = hgsfzero(f,x1,options)
%***********************************************************************************************************
%* HGS 1.3
%* By Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Usage:
%       [ xc,yc ] = HGSFZERO( f,x1,options )
%
% Inputs:
%   - f    -> function to solve.
%   - x1 -> Beginning of the search interval
%   - options -> Structure containing:
%       -> x2: End of the search interval
%       -> fchange: 
%       -> epsx: Solution tolerance
%       -> epsy: Function tolerance
%       -> maxite: Maximum number of iterations
%       -> info: Solver information
%
% Output:
%   - xc -> solution of the function.
%   - yc -> value of the function at the solution point.
%   - flag -> exit flag:
%       -> 1: Solution converged
%       -> -1: Solution not converged
%
% See also FSOLVE, FZERO
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% Get values from structure
x2 = options.x2;
fchange = options.fchange;
epsx = options.epsx;
epsy = options.epsy;
maxite = options.maxite;
info = options.info;

% Compute values of the function
y1=f(x1);
y2=f(x2);

if y1*y2>0 % uh ? no sign change, we can't proceed safely 
    flag=-1;
    xc=x1;
    yc=y1;
    return
end

for i=1:maxite
    if x2-x1 > fchange % fchange should be about 1
        xc=(x1+x2)/2;
        fc=1;
    else
        xc=x1-y1*(x2-x1)/(y2-y1);
        fc=0;
    end
    
    yc=f(xc);
    if ~isempty(info)
        fprintf('myfzero i=%d fc=%d x1=%7.2e y1=%7.2e x2=%7.2e y2=%7.2e xc=%7.2e dx=%7.2e |yc|=%7.2e \n',...
            i,fc,x1,y1,x2,y2,xc,x2-x1,abs(yc) );
    end
    if abs(yc)<epsy && x2-x1<epsx
        break;
    end
    if yc*y1>0
        y1=yc;
        x1=xc;
        %flag=1;
    else
        y2=yc;
        x2=xc;
        %flag=2;
    end

end

% Exit flag
flag = 1;
% Error flag
if i == maxite, flag = -1; end

end