function [obj] = linError(obj,R,normRatio,evolRatio)
% linError - computes the linearization error
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    obj - nonlinear system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      03-January-2008 
% Last update:  15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%determine maximum allowed enlargement of reachable set
%compute travelled vector of zonotope center
x0=center(R);
eAt=obj.taylor.eAt;
A=obj.A;
f0=obj.linError.f0;

trans=(eAt-I)*x0+inv(A)*(eAt-I)*f0;


maxHull=obj.maxDer*options.timeStep;
maxHull=interval(infimum(maxHull),supremum(maxHull));

%compute interval of reachable set
IH=interval(R);

%compute intervals of total reachable set
IHtotal=IH+maxHull;
totalInt=interval(IHtotal);
%compute intervals of input
IHinput=interval(options.U);
inputInt=interval(IHinput);

%translate intervals by linearization point
deltaX=obj.linError.p.x;
deltaU=obj.linError.p.u;

totalInt=totalInt-deltaX;
inputInt=inputInt-deltaU;

%compute linearization error by passing the intervals to the Lagrange
%remainder mFile
error=lagrangeRemainder(totalInt,inputInt);
IHerror=interval(infimum(error),supremum(error));
obj.linError.error=zonotope(IHerror); %add interval to constant value due to the linearization

%------------- END OF CODE --------------