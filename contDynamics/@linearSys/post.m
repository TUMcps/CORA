function [Rnext,options] = post(linsys,~,params,options)
% post - computes the reachable continuous set for one time step
%
% Syntax:
%    [Rnext,options] = post(linsys,~,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% References:
%    [1] A. Girard, C. Le Guernic, and O. Maler, "Efficient computation of 
%        reachable sets of linear time-invariant systems with inputs" in 
%        Hybrid Systems: Computation and Control, ser. LNCS 3927. Springer, 
%        2006, pp. 257--271.
%    [2] A. Girard, "Reachability of uncertain linear systems using 
%        zonotopes" in Hybrid Systems: Computation and Control, 
%        ser. LNCS 3414. Springer, 2005, pp. 291–305.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-May-2007 
% Last update:   27-April-2009
%                29-June-2009
%                08-August-2011
%                25-July-2016 (intervalhull replaced by interval)
%                07-July-2018 (wrapping-free approach is no longer the default option)
%                27-July-2018 (changing input updated)
%                18-September-2018 (changing input again updated)
%                07-November-2018 (arguments of function updated)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%load data
eAt=linsys.taylor.eAt;

% check whether standard or wrapping-free approach should be used
if strcmp(options.linAlg,'wrapping-free')
    % option 1 (wrapping free) [1, Alg. 2]
    if isfield(options,'uTransVec') % check whether input vector changes
        [Rhom,Rhom_tp,Rtrans,inputCorr] = aux_inputInducedUpdates(linsys,params,options);
    else
        Rtrans = options.Rtrans;
        inputCorr = 0;
        Rhom=eAt*options.Rhom + center(Rtrans);
        if options.compTimePoint
            Rhom_tp=eAt*options.Rhom_tp + center(Rtrans);
        end
    end
    Raux=eAt*options.Raux;
    Rpar=options.Rpar + interval(Raux) + interval(Rtrans) + (-center(Rtrans));

elseif strcmp(options.linAlg,'standard')
    % option 2 (not wrapping-free) [2, Alg. 1]
    if isfield(options,'uTransVec') % check whether input vector changes
        [Rhom,Rhom_tp,~,inputCorr] = aux_inputInducedUpdates(linsys,params,options);
    else
        Rtrans = options.Rtrans;
        inputCorr = 0;
        Rhom=eAt*options.Rhom + Rtrans;

        % computation of time-point solution desired?
        if options.compTimePoint
            Rhom_tp=eAt*options.Rhom_tp + Rtrans;
        end
    end
    % propagate particular solution
    Raux = eAt*options.Raux;

    % reduction operation
    Rpar = options.Rpar + Raux;
    if ~isnumeric(options.Rpar)
        Rpar = reduce(Rpar,options.reductionTechnique,options.zonotopeOrder);
    end
    Rhom = reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);

    % computation of time-point solution desired?
    if options.compTimePoint
        Rhom_tp=reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
    end
else
    throw(CORAerror('CORA:wrongFieldValue','options.alg',...
        {'standard','wrapping-free'}));
end

%save homogeneous and particulate solution to options struct
options.Rhom=Rhom;
if options.compTimePoint
    options.Rhom_tp=Rhom_tp;
end
options.Rpar=Rpar;
options.Raux=Raux;

%write results to reachable set struct Rnext
if isa(Rhom,'polytope')
    Rnext.ti=Rhom+polytope(Rpar)+polytope(inputCorr);
    if options.compTimePoint
        Rnext.tp=Rhom_tp+polytope(Rpar);
    else
        Rnext.tp = [];
    end
else
    Rnext.ti=Rhom+zonotope(Rpar)+inputCorr;
    if options.compTimePoint
        Rnext.tp=Rhom_tp+zonotope(Rpar);
    else
        Rnext.tp = [];
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [Rhom,Rhom_tp,Rtrans,inputCorr] = aux_inputInducedUpdates(linsys,params,options)
% recalculates set if options.uTrans has changed

eAt = linsys.taylor.eAt;
eAtInt = linsys.taylor.eAtInt;
F = linsys.taylor.F;
inputF = linsys.taylor.inputF;

% solution due to constant inputs
vTrans = zonotope(linsys.B*params.uTrans + center(linsys.E*params.W) + linsys.c);

Rtrans = eAtInt*vTrans;
% effect should only be considered once in a single time interval:
inputCorr = inputF*vTrans;

% homogeneous time-point solution
Rinit = options.Rhom_tp;
if strcmp(options.linAlg,'wrapping-free')
    Rhom_tp = eAt*Rinit + center(Rtrans);
elseif any(strcmp(options.linAlg,{'standard','fromStart'}))
    Rhom_tp = eAt*Rinit + Rtrans;
end

% homogeneous time-interval solution
% note: inputCorr considered later
if isa(Rinit,'polyZonotope') 
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit);
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1};
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit;
end

end

% ------------------------------ END OF CODE ------------------------------
