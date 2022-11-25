function Rerror = errorSolution(obj,options,Vdyn,Vstat)
% errorSolution - computes the solution due to the linearization error
%
% Syntax:  
%    Rerror = errorSolution(obj,options,Vdyn,Vstat)
%
% Inputs:
%    obj - linearized system
%    options - options struct (for nonlinear system)
%    Vdyn - set of admissible errors (dynamic)
%    Vstat - set of admissible errors (static) (optional)
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
%               04-May-2020 (MW, unification with errorSolutionQuad)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 4 || isempty(Vstat)
    errorStat = 0;
else
    % including static error
    errorStat = obj.taylor.eAtInt * Vstat;
end

%load data from object/options structure
Apower=obj.taylor.powers;
E=obj.taylor.error;
taylorTerms=options.taylorTerms;
r=options.timeStep;
factors = options.factor;

%initialize Asum
Asum=r*Vdyn;

for i=1:taylorTerms
    %compute powers
    ApowerV=factors(i+1)*Apower{i}*Vdyn;
    %compute sums
    Asum=Asum+ApowerV;
end

%get error due to finite Taylor series
if isa(Vdyn,'zonoBundle')
    F=E*Vdyn.Z{1}*r;
else
    F=E*Vdyn*r;
end

%Compute error solution (dyn. + stat.)
Rerror = Asum+F + errorStat;


%------------- END OF CODE --------------