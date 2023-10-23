function I = interval(Z)
% interval - overapproximates a zonotope by an interval
%
% Syntax:
%    I = interval(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    I - interval object
%
% Example: 
%    Z = zonotope([-1;1],[3 2 -1; 2 1 -1]);
%    I = interval(Z);
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(I,[1,2],'r');
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       14-September-2006 
% Last update:   22-March-2007
%                11-November-2010
%                04-May-2011
%                22-July-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract center
c = Z.c;

% determine lower and upper bounds in each dimension
delta = sum(abs(Z.G),2);
leftLimit = c - delta;
rightLimit = c + delta;

% instantiate interval
I = interval(leftLimit,rightLimit);

% ------------------------------ END OF CODE ------------------------------
