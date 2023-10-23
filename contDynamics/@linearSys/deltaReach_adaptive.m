function [Rdelta,options] = deltaReach_adaptive(obj,Rinit,options)
% deltaReach_adaptive - computes the reachable continuous set of
% the difference to  the initial state for all initial states
% implemented only in adaptive nonlinear analysis: no reduction,
% less options for chosen set representation
%
% Syntax:
%    [Rdelta,options] = deltaReach_adaptive(obj,Rinit,options)
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       15-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute delta reachable set
%load data from object structure
eAt=obj.taylor.eAt;
F=obj.taylor.F;
inputCorr=obj.taylor.inputCorr;
Rtrans=obj.taylor.Rtrans;

%first time step homogeneous solution
dim = length(F);
Rhom_tp_delta = (eAt - eye(dim))*Rinit + Rtrans;

O = polyZonotope(zeros(dim,1),[],[],[]);
% preliminary Rdelta (without RV)
Rdelta = enclose(O,Rhom_tp_delta) + F*zonotope(Rinit) + inputCorr;

% no reduction

%original computation
if options.isRV
    Rdelta = Rdelta + obj.taylor.RV;
end


% ------------------------------ END OF CODE ------------------------------
