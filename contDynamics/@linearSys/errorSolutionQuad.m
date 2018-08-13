function [Rerror] = errorSolutionQuad(obj,Vstat,Vdyn,options)
% errorSolutionQuad - computes the solution due to the linearization error
% in a specialized version when the input is a quadZonotope
%
% Syntax:  
%    [Rerror] = errorSolutionQuad(obj,V,options)
%
% Inputs:
%    obj - linear system object
%    V - quadZonotope object
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
% Written:      05-September-2012
% Last update:  18-September-2012
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%load data from object/options structure
Apower=obj.taylor.powers;
E=obj.taylor.error;
eAtInt=obj.taylor.eAtInt;
taylorTerms=options.taylorTerms;
r=options.timeStep;
dim=dimension(obj);

%initialize Asum
Asum=eye(dim)*r*Vdyn;

for i=1:taylorTerms
    %compute powers
    ApowerV{i}=options.factor(i+1)*Apower{i}*Vdyn;   
    %compute sums
    Asum=Asum+ApowerV{i};
end

%get error due to finite Taylor series
F=E*Vdyn*r;

%Compute error solution
Rerror=eAtInt*Vstat + Asum + F;

%Convert to zonotope if inputSol is an interval hull
if strcmp('interval',class(Rerror))
    Rerror=zonotope(interval(infimum(Rerror),supremum(Rerror)));
end

%------------- END OF CODE --------------