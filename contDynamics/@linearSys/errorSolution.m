function [Rerror] = errorSolution(obj,V,options)
% errorSolution - computes the solution due to the linearization error
%
% Syntax:  
%    [Rerror] = errorSolution(linSys,options)
%
% Inputs:
%    obj - linear system object
%    nonlinSys - nonlinear system object
%    options - options struct
%
% Outputs:
%    Rerror - reachable set due to the linearization error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      30-October-2007 
% Last update:  22-January-2008
%               18-May-2011
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object/options structure
Apower=obj.taylor.powers;
E=obj.taylor.error;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);

%initialize Asum
Asum=eye(dim)*r*V;

for i=1:taylorTerms
    %compute powers
    ApowerV{i}=r^(i+1)/factorial(i+1)*Apower{i}*V;   
    %compute sums
    Asum=Asum+ApowerV{i};
end

%get error due to finite Taylor series
if isa(V,'zonotopeBundle')
    F=E*V.Z{1}*r;
else
    F=E*V*r;
end

%Compute error solution
Rerror=Asum+F;

%Convert to zonotope if inputSol is an interval hull
if strcmp('interval',class(Rerror))
    Rerror=zonotope(interval(infimum(Rerror),supremum(Rerror)));
end

%------------- END OF CODE --------------