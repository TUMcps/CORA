function [Rfirst, options] = rA_initReach(obj, options)
% rA_initReach - init step: R([0, Delta t_0])
%
% Syntax:  
%    [Rfirst, options] = rA_initReach(obj, options)
%
% Inputs:
%    obj      - linearSys object
%    options  - options struct, containing
%                    .R0: initial set
%                    .timeStep: time step size for current iteration
%                    .taylorTerms: taylor Terms for current iteration
%
% Outputs:
%    Rfirst  - reachable set for time interval (.ti) and time point (.tp)
%    options - options struct, containing
%                    .Rhom: homogeneous time interval solution
%                    .Rhom_tp: homogeneous time point solution
%                    .isInhom: if inhomogenuity present
%                    .Rtrans: inhomogeneous solution due to uTrans
%                    .Raux: inhomogeneous solution due to U
%                    .Rinhom: full inhomogeneous solution
%                    .P: identity matrix (used for propagation of Rtrans)
%                    .Q: propagation matrix using .timeStep (for Raux)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% preliminary computations ------------------------------------------------
% compute exponential matrix:
% returns: obj.taylor.powers|error
% obj = exponential(obj,options);

% compute time interval error (tie):
% returns: obj.taylor.F
% obj = tie(obj,options);

% ^ both above already calculated in init_TimeStepTaylorTerms

if options.isInput
    options.factor = options.savedSets{1}.factor;
end
% compute reachable set due to input:
% returns: obj.taylor.V|RV|Rtrans|inputF|inputCorr|eAtInt
obj = inputSolution(obj,options);

% save time step using for above calculations in same struct
obj.taylor.timeStep = options.timeStep;

% compute reachable set of first time interval
eAt = expm(obj.A*options.timeStep);

% initialize propagation matrices for analysis
options.P = speye(obj.dim);
options.Q = eAt;

% save current propagation matrix in same struct
obj.taylor.eAt = eAt;
% -------------------------------------------------------------------------


% read necessary sets from preliminary computations -----------------------
F = obj.taylor.F;
Rinit = options.R0;
if any(any(options.U.Z))
    gensRV = generators(obj.taylor.RV);
    RV = zonotope([zeros(obj.dim,1), gensRV(:,any(gensRV,1))]);
else
    RV = zonotope(zeros(obj.dim,1));
end
if any(options.uTrans)
    if iscell(obj.taylor.Rtrans)
        Rtrans = obj.taylor.Rtrans{1};
    else
        gensRtrans = generators(obj.taylor.Rtrans);
        Rtrans = zonotope([center(obj.taylor.Rtrans), gensRtrans(:,any(gensRtrans,1))]);
    end
else
    Rtrans = zonotope(zeros(obj.dim,1));
end
inputCorr = obj.taylor.inputCorr;
% -------------------------------------------------------------------------


% first time step homogeneous solution ------------------------------------
Rhom_tp = eAt*Rinit + Rtrans;
if isa(Rinit,'quadZonotope')
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit) + inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1} + inputCorr;
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit + inputCorr;
end
% note: no reduction
% -------------------------------------------------------------------------


% save homogeneous and inhomogeneous solution -----------------------------
options.Rhom    = Rhom;
options.Rhom_tp = Rhom_tp;

options.Raux    = RV;
options.Rtrans  = Rtrans;
options.lastRtrans = options.Rtrans;

options.Rinhom = RV;
if options.isInput
    % generators of Rinhom
    options.RinhomG = generators(options.Rinhom);
    options.RinhomGlength = vecnorm(options.RinhomG,2);
    options.RinhomRed = 0;
    % propInhom_new ---
%     options.RinhomGorder = vecnorm(options.RinhomG,1) - vecnorm(options.RinhomG,Inf);
    % ---
    % container for generators converted to interval (always diag matrix)
    options.RinhomConvInt = zeros(obj.dim);
end
% -------------------------------------------------------------------------


% total solution ----------------------------------------------------------
if isa(Rinit,'mptPolytope')
    % convert zonotopes to polytopes
    Radd      = mptPolytope(RV);
    Rfirst.ti = Rhom + Radd;
    Rfirst.tp = Rhom_tp + Radd;
else
    % original computation
    Rfirst.ti = Rhom + RV;
    Rfirst.tp = Rhom_tp + RV;
end
% -------------------------------------------------------------------------



end


%------------- END OF CODE --------------

