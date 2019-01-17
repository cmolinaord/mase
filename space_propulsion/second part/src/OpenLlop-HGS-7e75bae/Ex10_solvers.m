%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 10: Solving an arbitrary equation with different solvers, calling
% hsgfsolve

clear
clc
close all

eq=@(x) sin(x)-exp(-x);


% solve with fsolve, starting iterations with 0.7; display progress and set
% x tolerance = 1e-10
[xc1,yc1]=hgssolve(eq,0.7,'fsolve',optimset('Display','iter','TolX',1e-6))

% the same, but with default options
[xc2,yc2]=hgssolve(eq,0.7,'fsolve',[])


% with fzero
% fzero begins searching an interval around 0.7 where sign of eq changes
% it is slower but more robust
[xc3,yc3]=hgssolve(eq,0.7,'fzero',optimset('Display','iter'))


% with hgsfzero
% this in-house code looks for the solution in the interval x1,x2 (where x1
% is the first argument and x2 is specified by the options). The function
% will never be evaluated outside [x1,x2]. If there is no sign change,
% hgsfzero refuses to begin. When the search interval is less than fchange,
% the next iteration point is choosed using the secant method.
% use hgsfzero as the last option, it is the most robust, but the slowest

% it will work with:
% eq(0.1)<0
% eq(0.6)>0

% info=0 -> don't display progress 
[xc4,yc4]=hgssolve(eq,0.1,'hgsfzero',struct('x2',0.6,'fchange',1e-2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',1));

% But it will not work with 0.1 & 0.2
[xc4,yc4]=hgssolve(eq,0.1,'hgsfzero',struct('x2',0.2,'fchange',1e-2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',1));


