function [trueError,VerrorDyn] = abstrerr_lin(obj,options,R)
% abstrerr_lin - computes the abstraction error for linearization approach;
%    to enter, options.alg = 'lin' and options.tensorOrder = 2|3
%
% Syntax:  
%    [trueError,VerrorDyn] = abstrerr_lin(obj,options,R)
%
% Inputs:
%    obj - nonlinearSys or nonlinParamSys object
%    options - options struct
%    R - reachable set (time-interval solution from linearized system
%           + estimated set of abstraction errors)
%
% Outputs:
%    trueError - abstraction error (interval)
%    Verrordyn - abstraction error (zonotope)
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   [1] M. Althoff et al. "Reachability Analysis of Nonlinear Systems with 
%       Uncertain Parameters using Conservative Linearization"
%
% See also: linReach
%
% Former files: linError.m, linError_higherOrder.m

% Author:        Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:       21-April-2020
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% compute interval of reachable set
IHx = interval(R);
% compute intervals of total reachable set
totalInt_x = IHx + obj.linError.p.x;

% compute intervals of input
IHu = interval(options.U);
% translate intervals by linearization point
totalInt_u = IHu + obj.linError.p.u;


if options.tensorOrder == 2
    
    % assign correct hessian (using interval arithmetic)
    obj = setHessian(obj,'int');

    % evaluate the hessian matrix with the selected range-bounding technique
    if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
       ~strcmp(options.lagrangeRem.method,'interval')

        % create taylor models or zoo-objects
        [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

        % evaluate the Lagrange remainder 
        if isa(obj,'nonlinParamSys')
            H = obj.hessian(objX,objU,options.paramInt);
        else
            H = obj.hessian(objX,objU);
        end
    else
        if isa(obj,'nonlinParamSys')
            H = obj.hessian(totalInt_x,totalInt_u,options.paramInt);
        else
            H = obj.hessian(totalInt_x,totalInt_u);
        end
    end

    if true
        % standard method using interval arithmetic

        % obtain maximum absolute values within IHx, IHu
        dx = max(abs(infimum(IHx)),abs(supremum(IHx)));
        du = max(abs(infimum(IHu)),abs(supremum(IHu)));
    
        % calculate the Lagrange remainder (second-order error)
        % ...acc. to Proposition 1 in [1]
        errorLagr = zeros(length(H),1);
        dz = [dx;du];
    
        for i = 1:length(H)
            H_ = abs(H{i});
            H_ = max(infimum(H_),supremum(H_));
            errorLagr(i) = 0.5 * dz' * H_ * dz;
        end
        
        trueError = errorLagr;
        VerrorDyn = zonotope([0*trueError,diag(trueError)]);

    else
        % no interval arithmetic (?)
        % attention: requires options.errorOrder|intermediateOrder!
%         options.errorOrder = 20;
%         options.intermediateOrder = 10;

        %compute zonotope of state and input
        Rred = reduce(zonotope(R),options.reductionTechnique,options.errorOrder);
        Z = cartProd(Rred,options.U);

        %obtain combined interval z and absolute values
        dz = [IHx; IHu];
        dz_abs = max(abs(infimum(dz)),abs(supremum(dz)));
        
        %separate evaluation
        H_mid = cell(obj.dim,1);
        H_rad = cell(obj.dim,1);
        for i=1:obj.dim
            H_mid{i} = sparse(center(H{i}));
            H_rad{i} = sparse(rad(H{i}));
        end
        error_mid = 0.5*quadMap(Z,H_mid);
        
        %interval evaluation
        error_rad = zeros(obj.dim,1);
        for i=1:obj.dim
            error_rad(i,1) = 0.5*dz_abs'*H_rad{i}*dz_abs;
        end
        
        %combine results
        VerrorDyn = error_mid + zonotope(interval(-error_rad,error_rad));
        VerrorDyn = reduce(VerrorDyn,options.reductionTechnique,options.intermediateOrder);
        
        % interval for containment check of subsequent Lagrange remainders
        trueError = supremum(abs(interval(VerrorDyn)));

    end    

elseif options.tensorOrder == 3

    % set handles to correct files
    obj = setHessian(obj,'standard');
    obj = setThirdOrderTensor(obj,'int');
    
    % obtain intervals and combined interval z
    dz = [IHx; IHu];

    % reduce zonotope
    Rred = reduce(R,options.reductionTechnique,options.errorOrder);

    % combined zonotope (states + input)
    Z = cartProd(Rred,options.U);

    % calculate hessian matrix
    if isa(obj,'nonlinParamSys')
        H = obj.hessian(obj.linError.p.x,obj.linError.p.u,options.paramInt);
    else
        H = obj.hessian(obj.linError.p.x,obj.linError.p.u);
    end

    % evaluate third-order tensor 
    if isfield(options,'lagrangeRem') && isfield(options.lagrangeRem,'method') && ...
        ~strcmp(options.lagrangeRem.method,'interval')

        % create taylor models or zoo objects
        [objX,objU] = initRangeBoundingObjects(totalInt_x,totalInt_u,options);

        if isa(obj,'nonlinParamSys')
            [T,ind] = obj.thirdOrderTensor(objX, objU, options.paramInt);
        else
            [T,ind] = obj.thirdOrderTensor(objX, objU);
        end

    else
        if isa(obj,'nonlinParamSys')
            [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u, options.paramInt);
        else
            [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u);
        end
    end

    % second-order error
    errorSec = 0.5 * quadMap(Z,H);

    % calculate the Lagrange remainder term (third-order error)
    if isempty(cell2mat(ind))
        % empty zonotope if all entries in T are empty
        errorLagr = zonotope(zeros(obj.dim,1));
    else
        % skip tensors with all-zero entries using ind from tensor creation
        errorLagr = interval(zeros(obj.dim,1),zeros(obj.dim,1));
        for i=1:length(ind)
            error_sum = interval(0,0);
            for j=1:length(ind{i})
                error_tmp = dz.'*T{i,j}*dz;
                error_sum = error_sum + error_tmp * dz(j);
            end
            errorLagr(i,1) = 1/6*error_sum;
        end
        errorLagr = zonotope(errorLagr);
    end

    % overall linearization error
    VerrorDyn = errorSec + errorLagr;
    VerrorDyn = reduce(VerrorDyn,...
        options.reductionTechnique,options.intermediateOrder);

    trueError = supremum(abs(interval(VerrorDyn)));
    
else
    throw(CORAerror('CORA:notSupported',...
        "No abstraction error computation for chosen tensor order!"));
end

%------------- END OF CODE --------------