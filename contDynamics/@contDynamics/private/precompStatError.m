function [H,Zdelta,errorStat,T,ind3,Zdelta3] = precompStatError(obj,Rdelta,options)
% precompStatError - precompute the second order static error along with 
%                    hessian matrix
%
% Syntax:  
%    [H,Zdelta,errorStat,T,ind3] = precompStatError(obj,Rdelta,options)
%
% Inputs:
%    obj - nonlinear system object
%    Rdelta - shifted reachable set at the beginning of the time step
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
% See also: 

% Author:       Niklas Kochdumper
% Written:      27-July-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

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
Ustat = zonotope(zeros(dim(options.U),1));
Z = cartProd(Rred,Ustat);
Zdelta = cartProd(Rdelta,Ustat);

% calculate the hessian tensor
if isa(obj,'nonlinParamSys')
    H = obj.hessian(obj.linError.p.x, obj.linError.p.u,options.paramInt);
else
    H = obj.hessian(obj.linError.p.x, obj.linError.p.u);
end

% calcualte the quadratic map == static second order error
errorSecOrdStat = 0.5*quadMap(Z, H);

% third-order error
if options.tensorOrder >= 4
   
    % reduce the order of the reachable set to speed-up the computations 
    % for cubic multiplication
    if isfield(options,'errorOrder3')
        
       Rred = reduce(Rred,options.reductionTechnique,options.errorOrder3);
       Rdelta = reduce(Rdelta,options.reductionTechnique,options.errorOrder3);
       
       Z = cartProd(Rred,options.U);
       Zdelta3 = cartProd(Rdelta,options.U);
       
    else
       Zdelta3 = Zdelta;
    end
    
    % calculate the third-order tensor
    if isa(obj,'nonlinParamSys')
        [T,ind3] = obj.thirdOrderTensor(obj.linError.p.x, ...
                                        obj.linError.p.u, options.paramInt);
    else
        [T,ind3] = obj.thirdOrderTensor(obj.linError.p.x, obj.linError.p.u);
    end
    
    % calculate the cubic map == static third-order error
    errorThirdOrdStat = 1/6 * cubMap(Z,T,ind3);
    
    % calculate the overall static linearization error
    if isa(Z,'polyZonotope')
        errorStat = exactPlus(errorSecOrdStat,errorThirdOrdStat);
    else
        errorStat = errorSecOrdStat + errorThirdOrdStat;
    end
    
else
    
    errorStat = errorSecOrdStat;
end

% reduce the complexity of the set of static errors
errorStat = reduce(errorStat,options.reductionTechnique,options.intermediateOrder);

%------------- END OF CODE -------------