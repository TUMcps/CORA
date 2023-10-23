function [Rdelta,options] = deltaReach(obj,Rinit,options)
% deltaReach - computes the reachable continuous set of the difference to 
% the initial state for all initial states
%
% Syntax:
%    [Rdelta,options] = deltaReach(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
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

% Authors:       Matthias Althoff
% Written:       18-September-2012
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%compute delta reachable set
%load data from object structure
eAt=obj.taylor.eAt;
F=obj.taylor.F;
RV=obj.taylor.RV;
inputCorr=obj.taylor.inputCorr;
Rtrans=obj.taylor.Rtrans;


%first time step homogeneous solution
dim = length(F);
Rhom_tp_delta = (eAt - eye(dim))*Rinit + Rtrans;

if isa(Rinit,'zonotope')
    %original computation
    O = zonotope(zeros(dim,1));
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit+inputCorr;
elseif isa(Rinit,'polyZonotope') || isa(Rinit,'conPolyZono')
    O = zeros(dim)*Rhom_tp_delta;
    Rhom=enclose(O,Rhom_tp_delta)+F*zonotope(Rinit)+inputCorr;
elseif isa(Rinit,'zonoBundle')
    O = zonoBundle({zonotope(zeros(dim,1))});
    Rhom=enclose(O,Rhom_tp_delta)+F*Rinit.Z{1}+inputCorr;
end

%reduce zonotope
Rhom=reduce(Rhom,options.reductionTechnique,options.intermediateOrder);
if ~isnumeric(RV)
    RV=reduce(RV,options.reductionTechnique,options.intermediateOrder);
end


%total solution
if isa(Rinit,'polytope')
    %convert zonotopes to polytopes
    Radd=polytope(RV);
    Rdelta=Rhom+Radd;
else
    %original computation
    Rdelta=Rhom+RV;
end


% ------------------------------ END OF CODE ------------------------------
