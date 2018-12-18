function [Rdelta,options] = deltaReach(obj,Rinit,options)
% deltaReach - computes the reachable continuous set of the difference to 
% the initial state for all initial states
%
% Syntax:  
%    [obj,Rfirst,options] = deltaReach(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%    Rdelta - difference reachable set 
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written:      18-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%compute delta reachable set
%load data from object structure
eAt=obj.taylor.eAt;
F=obj.taylor.F;
RV=obj.taylor.RV;
inputCorr=obj.taylor.inputCorr;
Rtrans=obj.taylor.Rtrans;


%first time step homogeneous solution
dim = length(F);
Rhom_tp_delta=(eAt - eye(dim))*Rinit + Rtrans;
%change in computation for polytope examples
if isa(Rinit,'mptPolytope')
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit;
elseif isa(Rinit,'quadZonotope')
    O = quadZonotope(zeros(dim,1));
    Rhom=enclose(O,Rhom_tp_delta)+F*zonotope(Rinit)+inputCorr;
elseif isa(Rinit,'zonotopeBundle')
    O = zonotopeBundle({zonotope(zeros(dim,1))});
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit.Z{1}+inputCorr;
else
    %original computation
    O = zonotope(zeros(dim,1));
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit+inputCorr;
end

%reduce zonotope
Rhom=reduce(Rhom,options.reductionTechnique,options.intermediateOrder);
RV=reduce(RV,options.reductionTechnique,options.intermediateOrder);


%total solution
if isa(Rinit,'mptPolytope')
    %convert zonotopes to polytopes
    Radd=mptPolytope(RV);
    Rdelta=Rhom+Radd;
else
    %original computation
    Rdelta=Rhom+RV;
end


%------------- END OF CODE --------------