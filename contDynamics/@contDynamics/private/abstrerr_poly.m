function [trueError, VerrorDyn, VerrorStat] = ...
    abstrerr_poly(obj, options, Rall, Rdiff, ...
    H, Zdelta, VerrorStat, T, ind3, Zdelta3)
% abstrerr_poly - computes the abstraction error for the polynomialization
%    approach introduced in [1]
%
% Syntax:  
%    [errorStat, errorDyn, errorInt] = 
%           linError_thirdOrder(obj, options, Rall, Rdiff, H, Zdelta, 
%                               errorStat, T, ind3, Zdelta3)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    Rall - time-interval reachable set
%    Rdiff - difference between the reachable set at the beginning of the
%            time interval and the time-interval reachable set
%    (the following six input parameters come from precompStatError.m)
%    H - Hessian matrix
%    Zdelta - zonotope over-approximating the reachable set at the
%             beginning of the time step extended by the input set
%    errorStat - set of static linearization errors
%    T - third-order tensor
%    ind3 - indices of non-zero entries in the third-order tensor
%    Zdelta3 - set Zdelta reduced to the zonotope order for the evaluation
%              of the third-order tensor
%
% Outputs:
%    trueError - interval overapproximating the overall linearization error 
%    errorDyn - zonotope overapproximating the dynamic linearization error
%    errorStat - zonotope overapproximating the static linearization error
%
% Other m-files required: precompStatError.m
% Subfunctions: none
% MAT-files required: none
%
% See also: linError_higherOrder.m
%
% References: 
%   [1] M. Althoff
%       "Reachability Analysis of Nonlinear Systems using
%           Conservative Polynomialization and Non-Convex Sets"

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  18-March-2016
%               25-July-2016 (intervalhull replaced by interval)
%               22-January-2018 (NK, fixed error for the sets)
%               08-February-2018 (NK, higher-order-tensors + clean-up)
%               21-April-2020 (simplification)
% Last revision:---

%------------- BEGIN CODE --------------

% compute interval of reachable set
dx = interval(Rall);
totalInt_x = dx + obj.linError.p.x;

% compute intervals of input
du = interval(options.U);
totalInt_u = du + obj.linError.p.u;

% obtain intervals and combined interval z
dz = [dx; du];

% compute zonotope of state and input
Rred_diff = reduce(zonotope(Rdiff),options.reductionTechnique,options.errorOrder);
Z_diff = cartProd(Rred_diff,options.U);

% second-order error
error_secondOrder_dyn = 0.5*(quadMap(Zdelta,Z_diff,H) ...
                             + quadMap(Z_diff,Zdelta,H) + quadMap(Z_diff,H));

% third-order error
if options.tensorOrder == 3
    
    % evaluate the third-order tensor
    if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
       ~strcmp(options.lagrangeRem.method,'interval')

        % create taylor models or zoo-objects
        [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

        % evaluate third order tensor 
        if isa(obj,'nonlinParamSys')
            [T,ind] = obj.thirdOrderTensor(objX, objU, options.paramInt);
        else
            [T,ind] = obj.thirdOrderTensor(objX, objU);
        end

    else
        if isa(obj,'nonlinParamSys')
            [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u,options.paramInt);
        else
            [T,ind] = obj.thirdOrderTensor(totalInt_x,totalInt_u);
        end
    end
    
    % calculate the Lagrange remainder term
    error_thirdOrder_dyn = interval(zeros(obj.dim,1),zeros(obj.dim,1));
    for i=1:length(ind)
        error_sum = interval(0,0);
        for j=1:length(ind{i})
            error_sum = error_sum + (dz.'*T{i,ind{i}(j)}*dz) * dz(ind{i}(j));
        end
        error_thirdOrder_dyn(i,1) = 1/6*error_sum;
    end

    error_thirdOrder_dyn = zonotope(error_thirdOrder_dyn);
    
    % no terms of order >= 4
    remainder = zonotope(zeros(obj.dim,1));
    
else
    % tensorOrder >= 4
    
    % reduce set Zdiff to the desired zonotope order to speed up the
    % computation of cubic multiplication
    if isfield(options,'errorOrder3')
         Rred_diff = reduce(Rred_diff,...
             options.reductionTechnique,options.errorOrder3);
         Z_diff3 = cartProd(Rred_diff,options.U);
    else
         Z_diff3 = Z_diff;
    end
    
    % third-order error
    error_thirdOrder_dyn = 1/6*(cubMap(Zdelta3,T,ind3) + ...
                                cubMap(Zdelta3,Zdelta3,Z_diff3,T,ind3) + ...
                                cubMap(Zdelta3,Z_diff3,Z_diff3,T,ind3) + ...
                                cubMap(Zdelta3,Z_diff3,Zdelta3,T,ind3) + ... 
                                cubMap(Z_diff3,Zdelta3,Z_diff3,T,ind3) + ...
                                cubMap(Z_diff3,Zdelta3,Zdelta3,T,ind3) + ...
                                cubMap(Z_diff3,Z_diff3,Zdelta3,T,ind3));
    
    % init higher-order error
    remainder = interval(zeros(obj.dim,1),zeros(obj.dim,1));

    % exact evaluation of intermediate taylor terms
    for i=4:options.tensorOrder-1
        handle = obj.tensors{i-3};
        remainder = remainder + handle(obj.linError.p.x,obj.linError.p.u,dx,du);
    end

    % lagrange remainder over-approximating the last taylor term
    handle = obj.tensors{options.tensorOrder-3};

    if isa(obj,'nonlinParamSys')
        remainder = remainder + handle(totalInt_x,totalInt_u,dx,du,options.paramInt);
    else
        remainder = remainder + handle(totalInt_x,totalInt_u,dx,du);
    end
    remainder = zonotope(remainder);
    
end

% combine results
VerrorDyn = error_secondOrder_dyn + error_thirdOrder_dyn + remainder;
VerrorDyn = reduce(VerrorDyn,options.reductionTechnique,options.intermediateOrder);

errorIHabs = abs(interval(VerrorDyn) + interval(VerrorStat));
trueError = supremum(errorIHabs);

%------------- END OF CODE --------------