function [Z_error, errorInt, errorInt_x, errorInt_y, R_y, options] = ...
            	abstractionError_adaptive(obj, options, R, Verror_y)
% abstractionError_adaptive - computes the abstraction error
% note: no Taylor Model or zoo integration
%
% Syntax:
%    [Verror, errorInt, errorInt_x, errorInt_y, R_y] = ...
%              abstractionError_adaptive(obj, options, R, Verror_y)
%
% Inputs:
%    obj - nonlinDASys object
%    options - options struct
%    R - reachable set
%    Verror_y - set of algebraic linearization error
%
% Outputs:
%    Z_error - zonotope over-approximating the linearization error
%    errorInt - interval over-approximating the linearization error
%    errorInt_x - interval over-approximating the linearization error (dynamic part)
%    errorInt_y - interval over-approximating the linearization error (constraint part)
%    R_y - reachable set of the algebraic part
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%    [1] M. Althoff, B. Krogh. "Reachability analysis of nonlinear
%        differential-algebraic systems", IEEE Transactions of
%        Automatic Control, 2013.
%
% See also: linReach, linError_*

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       30-August-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% read out variables
f0_con = obj.linError.f0_con;
D = obj.linError.D;
E = obj.linError.E;
F_inv = obj.linError.F_inv;

% intervals for reachable set and input
dx = interval(R);
du = interval(options.U);

% shift by linearization point
totalInt_x = dx + obj.linError.p.x;
totalInt_u = du + obj.linError.p.u;


% LIN; TENSOR 2 -----------------------------------------------------------
if strcmp(options.alg,'lin') && options.tensorOrder == 2
    % approach in [Sec. IV, 1]
    
    % set handle to correct file
    obj = setHessian(obj,'int');

    % compute set of algebraic variables
    F_inv = obj.linError.F_inv;
    R_y = -F_inv*(f0_con + D*R + E*options.U + Verror_y); % see [1, eq.(12)]

    % obtain interval of algebraic states and combined interval z
    dy = interval(R_y);
    dz = [dx; dy; du];

    % shift interval of algebraic states by linearization point
    totalInt_y = dy + obj.linError.p.y;
    
    % evaluate Hessian
    [Hf, Hg] = obj.hessian(totalInt_x, totalInt_y, totalInt_u);

    % error_x
    for i=1:length(Hf)
        Z_error_x(i,1) = 0.5*dz.'*Hf{i}*dz;
    end

    % error_y
    for i=1:length(Hg)
        Z_error_y(i,1) = 0.5*dz.'*Hg{i}*dz;
    end

    % final error
    Z_error_x = zonotope(interval(infimum(Z_error_x),supremum(Z_error_x)));
    Z_error_y = zonotope(interval(infimum(Z_error_y),supremum(Z_error_y)));
    Z_error = Z_error_x + obj.linError.CF_inv*Z_error_y;

    % update R_y, see [1, eq.(19)]
    R_y =  obj.linError.p.y + (-F_inv)*(f0_con + D*R + E*options.U + Z_error_y);

    % error intervals
    errorIHabs = abs(interval(Z_error));
    errorInt = supremum(errorIHabs);

    errorIHabs_y = abs(interval(Z_error_y));
    errorInt_y = supremum(errorIHabs_y);

    errorIHabs_x = abs(interval(Z_error_x));
    errorInt_x = supremum(errorIHabs_x);
    
    
% LIN; TENSOR 3 -----------------------------------------------------------
elseif strcmp(options.alg,'lin') && options.tensorOrder == 3
    
    % set handles to correct files
    obj = setHessian(obj,'standard');
    obj = setThirdOrderTensor(obj,'int');
    
    % compute set of algebraic variables, see [1, eq.(19)], split into two
    % parts acc. to [1, Prop. 2]
    R_y_cor = -F_inv*(f0_con + D*R); % x-y-correlated part
    R_y_add = -F_inv*(E*options.U + Verror_y); % uncorrelated part
    
    % obtain intervals and combined interval z
    dy = interval(R_y_cor + R_y_add);
    dz = [dx; dy; du];
    
    % shift interval of algebraic states by linearization point
    totalInt_y = dy + obj.linError.p.y;
    
    % evaluate Hessian (matrices) and third-order tensor (interval matrices)
    [Hf, Hg] = obj.hessian(obj.linError.p.x, obj.linError.p.y, obj.linError.p.u);
    [Tf, Tg] = obj.thirdOrderTensor(totalInt_x, totalInt_y, totalInt_u);
    
    % compute zonotope of state, constraint variables, and input
    Z_x = [R.c,R.G];
    Z_y_cor = [R_y_cor.c,R_y_cor.G];
    Z_y_add = [R_y_add.c,R_y_add.G];
    Z_0 = zeros(length(Z_x(:,1)), length(Z_y_add(1,:)));
    R_xy = zonotope([Z_x, Z_0; Z_y_cor, Z_y_add]);     % see [1, Prop. 2]
    R_xyu = cartProd(R_xy, options.U);
    % order reduction before application of quadMap operation
    if isfield(options,'gredIdx')
        if length(options.gredIdx.xyu) == options.i
            R_xyu = reduce(R_xyu,'idx',options.gredIdx.xyu{options.i});
        else
            [R_xyu,~,options.gredIdx.xyu{options.i,1}] = ...
                reduce(R_xyu,'adaptive',sqrt(options.redFactor));
        end
    else
        % safety branch (call from algorithm without gredIdx)
        R_xyu = reduce(R_xyu,'adaptive',sqrt(options.redFactor));
    end
    
    % second-order error for x and y (use quadMap_parallel instead of
    % quadMap to keep correlation... quadMap auto-removes 0-generators)
    error_x_secondOrder = 0.5*quadMap_parallel(R_xyu, Hf);
    error_y_secondOrder = 0.5*quadMap_parallel(R_xyu, Hg);
    
    % compute full second-order error
    % including correlation: exactplus (bug if quadMap used!)
    Z_err_x = [error_x_secondOrder.c,error_x_secondOrder.G];
    Z_err_x_add = obj.linError.CF_inv*error_y_secondOrder;
    Z_err_x_add = [Z_err_x_add.c,Z_err_x_add.G];
    error_secondOrder = zonotope(Z_err_x + Z_err_x_add);
    % remove 0-generators here... (remove zonotope-constructor above)
%     error_secondOrder = zonotope(error_secondOrder(:,any(error_secondOrder,1)));
    
    % third-order interval evaluation (dynamic part)
    for i=1:length(Tf(:,1))
        error_sum = interval(0,0);
        for j=1:length(Tf(1,:))
            error_tmp = dz'*Tf{i,j}*dz;
            error_sum = error_sum + error_tmp * dz(j);
        end
        error_x_thirdOrder(i,1) = 1/6*error_sum;
    end
    
    % third-order interval evaluation (algebraic part)
    for i=1:length(Tg(:,1))
        error_sum = interval(0,0);
        for j=1:length(Tg(1,:))
            error_tmp = dz'*Tg{i,j}*dz;
            error_sum = error_sum + error_tmp * dz(j);
        end
        error_y_thirdOrder(i,1) = 1/6*error_sum;
    end
    
    % convert to zonotopes
    error_thirdOrder_x_zono = zonotope(error_x_thirdOrder);
    error_thirdOrder_y_zono = zonotope(error_y_thirdOrder);
    
    % compute final error
    error_thirdOrder = error_thirdOrder_x_zono + ...
        obj.linError.CF_inv * error_thirdOrder_y_zono;
    Z_error = error_secondOrder + error_thirdOrder;
    
    % combine results
    Z_error_x = error_x_secondOrder + error_thirdOrder_x_zono;
    Z_error_y = error_y_secondOrder + error_thirdOrder_y_zono;
    
    % reduction
    if isfield(options,'gredIdx')
        if length(options.gredIdx.Z_error) == options.i
            Z_error = reduce(Z_error,'idx',options.gredIdx.Z_error{options.i});
            Z_error_y = reduce(Z_error_y,'idx',options.gredIdx.Z_error_y{options.i});
        else % also for options.run = 0 (initTensorOrder)
            [Z_error,~,options.gredIdx.Z_error{options.i,1}] = ...
                reduce(Z_error,'adaptive',10*options.redFactor);
            [Z_error_y,~,options.gredIdx.Z_error_y{options.i,1}] = ...
                reduce(Z_error_y,'adaptive',10*options.redFactor);
        end
    else
        % safety branch (call from algorithm without gredIdx)
        Z_error = reduce(Z_error,'adaptive',10*options.redFactor);
        Z_error_y = reduce(Z_error_y,'adaptive',10*options.redFactor);
    end
    
    % update R_y (time-interval)
    R_y = obj.linError.p.y + (-F_inv)*(f0_con + D*R + E*options.U + Z_error_y);
    
    % error intervals
    errorIHabs_x = abs(interval(Z_error_x));
    errorInt_x = supremum(errorIHabs_x);
    
    errorIHabs_y = abs(interval(Z_error_y));
    errorInt_y = supremum(errorIHabs_y);
    
    errorIHabs = abs(interval(Z_error));
    errorInt = supremum(errorIHabs);
    
    % OTHER COMBINATIONS OF ALG/TENSORORDER NOT IMPLEMENTED YET
else
    throw(CORAerror('CORA:notSupported','Given tensor order not supported.'));
end

end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkIfHessianConst(obj,options,H,totalInt_x,totalInt_u)
% check if hessian is constant --- only once executed! (very first step)
% returns: isHessianConst, hessianCheck, hessianConst

% assume true, check until false
options.isHessianConst = true;

% evaluate Hessians with larger set, check if result equal
H_test = obj.hessian(enlarge(totalInt_x,1.1),...
    enlarge(totalInt_u,1.1));
for i=1:length(H)
    if ~all(H{i} == H_test{i})
        options.isHessianConst = false; break
    end
end

% save only what is precisely needed
if options.isHessianConst
    for i=1:length(H)
        temp = abs(H{i});
        options.hessianConst{i} = ...
            max(infimum(temp),supremum(temp));
    end
end
options.hessianCheck = true;

end

% ------------------------------ END OF CODE ------------------------------
