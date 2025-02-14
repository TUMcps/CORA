function [H,Zdelta,errorStat,T,ind3,Zdelta3] = priv_precompStatError(sys,Rdelta,params,options)
% priv_precompStatError - precompute the second order static error along with 
%    Hessian matrix
%
% Syntax:
%    [H,Zdelta,errorStat,T,ind3] = priv_precompStatError(sys,Rdelta,params,options)
%
% Inputs:
%    sys - nonlinear system object
%    Rdelta - shifted reachable set at the beginning of the time step
%    params - model parameters
%    options - options struct
%
% Outputs:
%    H - hessian matrix
%    Zdelta - zonotope over-approximating the reachable set at the
%             beginning of the time step extended by the input set 
%    errorStat - static linearization error
%    T - third-order tensor
%    ind3 - indices at which the third-order tensor is not zero
%    Zdelta3 - set Zdelta reduced to the zonotope order for the evaluation
%              of the third-order tensor
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       27-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set handle to correct file
sys = setHessian(sys,'standard');

% initialize output arguments
T = [];
ind3 = [];
Zdelta3 = [];

% reduce the reachable set for the initial time point
Rred = reduce(Rdelta,options.reductionTechnique,options.errorOrder);

% compute a zonotope over-approximation of the reachable set at the
% initial time point
if isa(Rdelta,'zonotope')
   Rdelta = Rred;
else
   Rdelta = reduce(zonotope(Rdelta),options.reductionTechnique,options.errorOrder);
end

% extend the sets by the input sets
Ustat = zonotope(zeros(dim(params.U),1));
Z = cartProd(Rred,Ustat);
Zdelta = cartProd(Rdelta,Ustat);

% calculate the hessian tensor
if isa(sys,'nonlinParamSys')
    H = sys.hessian(sys.linError.p.x, sys.linError.p.u,params.paramInt);
else
    H = sys.hessian(sys.linError.p.x, sys.linError.p.u);
end

% calculate the quadratic map == static second order error
errorSecOrdStat = 0.5*quadMap(Z, H);

% third-order error
if options.tensorOrder >= 4
    
    % set handle to correct file
    sys = setThirdOrderTensor(sys,'standard');
   
    % reduce the order of the reachable set to speed-up the computations 
    % for cubic multiplication
    if isfield(options,'errorOrder3')
        
       Rred = reduce(Rred,options.reductionTechnique,options.errorOrder3);
       Rdelta = reduce(Rdelta,options.reductionTechnique,options.errorOrder3);
       
       Z = cartProd(Rred,params.U);
       Zdelta3 = cartProd(Rdelta,params.U);
       
    else
       Zdelta3 = Zdelta;
    end
    
    % calculate the third-order tensor
    if isa(sys,'nonlinParamSys')
        [T,ind3] = sys.thirdOrderTensor(sys.linError.p.x, ...
                                        sys.linError.p.u, params.paramInt);
    else
        [T,ind3] = sys.thirdOrderTensor(sys.linError.p.x, sys.linError.p.u);
    end
    
    % calculate the cubic map == static third-order error
    errorThirdOrdStat = 1/6 * cubMap(Z,T,ind3);
    
    % calculate the overall static linearization error
    if isa(Z,'polyZonotope') || isa(Z,'conPolyZono')
        errorStat = exactPlus(errorSecOrdStat,errorThirdOrdStat);
    else
        errorStat = errorSecOrdStat + errorThirdOrdStat;
    end
    
else
    
    errorStat = errorSecOrdStat;
end

% reduce the complexity of the set of static errors
errorStat = reduce(errorStat,options.reductionTechnique,options.intermediateOrder);

% ------------------------------ END OF CODE ------------------------------
