function [options, eAt] = rA_post(obj, options)
% rA_post - post step for adaptively parametrized reachability analysis
%
% Syntax: [options, eAt] = rA_post(obj, options)
%
% Inputs:
%    obj      - linearSys object
%    options  - options struct, containing
%                    .R0: initial set
%                    .timeStep: time step size for current iteration
%                    .taylorTerms: taylor Terms for current iteration
%
% Outputs:
%    options - options struct, containing
%                    .Rhom: homogeneous time interval solution
%                    .Rhom_tp: homogeneous time point solution
%                    .Rtrans: inhomogeneous solution due to uTrans
%                    .Raux: inhomogeneous solution due to U
%    eAt     - propagation matrix to end of current time step
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
% Last update:  31-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------

% preliminary computations ------------------------------------------------
% compute exponential matrix:
% returns: obj.taylor.powers|E
% obj = exponential(obj,options);

% compute time interval error (tie):
% returns: obj.taylor.F
% obj = tie(obj,options);

% ^ both above already calculated in rA_init_TimeStepTaylorTerms

% compute reachable set due to input:
% assert correct obj.taylor.error|powers
% returns: obj.taylor.V|RV|Rtrans|eAtInt
% ...also obj.taylor.inputF|inputCorr (already calculated in postTSTT...)
obj = inputSolution(obj,options);

% save time step using for above calculations in same struct
obj.taylor.timeStep = options.timeStep;

% compute reachable set of first time interval
eAt = expm(obj.A*options.timeStep);

% save current propagation matrix in same struct
obj.taylor.eAt = eAt;
% -------------------------------------------------------------------------


% read necessary sets from preliminary computations -----------------------
% assert obj.taylor.F|RV|Rtrans|inputCorr
F = obj.taylor.F;
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
        Rtrans = zonotope([center(obj.taylor.Rtrans), ...
            gensRtrans(:,any(gensRtrans,1))]);
    end
else
    Rtrans = zonotope(zeros(obj.dim,1));
end
inputCorr = obj.taylor.inputCorr;
Rinit = options.R0;
% -------------------------------------------------------------------------


% first time step homogeneous solution ------------------------------------
Rhom_tp = eAt*Rinit + Rtrans;
if isa(Rinit,'quadZonotope')
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit) + inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1} + inputCorr;
else
    enc = enclose(Rinit,Rhom_tp);
    bloat = F*Rinit;
    Rhom = enc + bloat + inputCorr;
end
% note: no reduction
% -------------------------------------------------------------------------

% save homogeneous and particulate solution -------------------------------
options.Rhom    = Rhom;
options.Rhom_tp = Rhom_tp;
options.Raux    = RV;
options.Rtrans  = Rtrans;
% -------------------------------------------------------------------------

end


%------------- END OF CODE --------------