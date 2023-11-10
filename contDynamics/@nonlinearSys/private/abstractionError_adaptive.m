function [VerrorDyn,VerrorStat,err,options] = abstractionError_adaptive(obj,options,R,Rdiff,...
    H,Zdelta,VerrorStat,T,ind3,Zdelta3)
% abstractionError_adaptive - computes the abstraction error
% note: no Taylor Model or zoo integration
%
% Syntax:
%    [VerrorDyn,VerrorStat,err,options] = abstractionError_adaptive(obj,options,R,Rdiff,...
%     H,Zdelta,VerrorStat,T,ind3,Zdelta3)
%
% Inputs:
%    obj - nonlinearSys object
%    options - options struct
%    R - time-interval solution of current step (incl. adm. abstr. err)
%    Rdiff - set of state differences [2,(6)]
%    Rdelta - initial set of step translated by linearization point
%    remaining inputs: same as outputs of precompStatError (only 'poly')
%
% Outputs:
%    VerrorDyn - set based on abstraction error
%    VerrorStat - set based on abstraction error
%    err - over-approximated abstraction error
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   [1] M. Althoff et al. "Reachability Analysis of Nonlinear Systems with 
%       Uncertain Parameters using Conservative Linearization"
%   [2] M. Althoff et al. "Reachability analysis of nonlinear systems using 
%       conservative polynomialization and non-convex sets"
%
% See also: linReach, linError_*, preCompStatError

% Authors:       Mark Wetzlinger
% Written:       14-January-2020
% Last update:   14-April-2020
%                20-November-2023 (MW, store generators selected for reduction)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute interval of reachable set and input (at origin)
IH_x = interval(R);
IH_u = interval(options.U);
% translate intervals by linearization point (at traj / lin. point)
totalInt_x = IH_x + obj.linError.p.x;
totalInt_u = IH_u + obj.linError.p.u;


% LIN; TENSOR 2 -----------------------------------------------------------
if strcmp(options.alg,'lin') && options.tensorOrder == 2
    % approach in [1]
    
    % assign correct hessian (using interval arithmetic)
    obj = setHessian(obj,'int');

    % obtain maximum absolute values within IH, IHinput
    IHinf = abs(infimum(IH_x));
    IHsup = abs(supremum(IH_x));
    dx = max(IHinf,IHsup);

    IHinputInf = abs(infimum(IH_u));
    IHinputSup = abs(supremum(IH_u));
    du = max(IHinputInf,IHinputSup);

    % compute an over-approximation of the Lagrange remainder [1,Prop.1]
    dz = [dx;du];
    
    % evaluate the hessian matrix with interval arithmetic
    try
        if ~options.isHessianConst
            H = obj.hessian(totalInt_x,totalInt_u);
            % very first step: check if Hessian is constant
            if ~options.hessianCheck
                options = aux_checkIfHessianConst(obj,options,H,totalInt_x,totalInt_u);
            end
        end
    catch ME
        if strcmp(ME.identifier,'interval:setoutofdomain')
            ME = MException('reach:setoutofdomain',...
                'Interval Arithmetic: Set out of domain.');
            throw(ME);
        else
            rethrow(ME);
        end
    end

    err = zeros(obj.dim,1);
    if ~options.isHessianConst
        for i = 1:obj.dim
            H_ = abs(H{i});
            H_ = max(infimum(H_),supremum(H_));
            err(i,1) = 0.5 * dz' * H_ * dz;
        end
    else
        % use saved H_ (the same in every step)
        for i = 1:obj.dim
            err(i,1) = 0.5 * dz' * options.hessianConst{i} * dz;
        end
    end

    VerrorDyn = zonotope(0*err,diag(err));
    VerrorStat = [];
    
% LIN; TENSOR 3 -----------------------------------------------------------
elseif strcmp(options.alg,'lin') && options.tensorOrder == 3

    % set handles to correct files
    obj = setHessian(obj,'standard');
    obj = setThirdOrderTensor(obj,'int');
    
    % calculate hessian matrix
    H = obj.hessian(obj.linError.p.x,obj.linError.p.u);
    
    % concatenate to interval z
    dz = [IH_x; IH_u];

    % reduce zonotope
    if isfield(options,'gredIdx')
        if length(options.gredIdx.Rred) == options.i
            % select the same generators for reduction as in the previous
            % iteration of the current time step to limit the influence of
            % the order reduction on the computation of the gain order
            if any(options.gredIdx.Rred{options.i} > size(generators(R),2))
                keyboard
            end
            Rred = reduce(R,'idx',options.gredIdx.Rred{options.i});
        else
            % reduce adaptively and store indices of generators that have
            % not been selected for reduction
            [Rred,~,options.gredIdx.Rred{options.i}] = ...
                reduce(R,'adaptive',sqrt(options.redFactor));
        end
    else
        % safety branch (call from algorithm without gredIdx)
        Rred = reduce(R,'adaptive',sqrt(options.redFactor));
    end
    
    % zonotope (states + input)
    Z = cartProd(Rred,options.U);
    
    % error second-order
    errorSec = 0.5 * quadMap(Z,H);
    
    % calculate third-order tensor
    try
        if ~options.thirdOrderTensorempty
            [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u);
        else
            ind = {};
        end
    catch ME
        if strcmp(ME.identifier,'interval:setoutofdomain')
            ME = MException('reach:setoutofdomain',...
                'Interval Arithmetic: Set out of domain.');
            throw(ME);
        else
            rethrow(ME);
        end
    end
    
    % calculate the Lagrange remainder term
    % skip tensors with all-zero entries using ind from tensor creation
    if ~isempty(cell2mat(ind))
        errorLagr = interval(zeros(obj.dim,1),zeros(obj.dim,1));
        for i=1:length(ind)
            error_sum = interval(0,0);
            for j=1:length(ind{i})
                error_tmp = dz.'*T{i,ind{i}(j)}*dz;
                error_sum = error_sum + error_tmp * dz(j);
            end
            errorLagr(i,1) = 1/6*error_sum;
        end
        errorLagr = zonotope(errorLagr);
    else
        errorLagr = 0;
        % third order tensor is empty, skip computation next time
        options.thirdOrderTensorempty = true;
    end
    
    % overall linearization error
    VerrorDyn = errorSec + errorLagr;
    if isfield(options,'gredIdx')
        if length(options.gredIdx.VerrorDyn) == options.i
            % select the same generators for reduction as in the previous
            % iteration of the current time step to limit the influence of
            % the order reduction on the computation of the gain order
            VerrorDyn = reduce(VerrorDyn,'idx',options.gredIdx.VerrorDyn{options.i});
        else
            % reduce adaptively and store indices of generators that have
            % not been selected for reduction
            [VerrorDyn,~,options.gredIdx.VerrorDyn{options.i}] = ...
                reduce(VerrorDyn,'adaptive',10*options.redFactor);
        end
    else
        % safety branch (call from algorithm without gredIdx)
        VerrorDyn = reduce(VerrorDyn,'adaptive',10*options.redFactor);
    end
    
    VerrorStat = [];
    
    err = abs(center(VerrorDyn)) + sum(abs(generators(VerrorDyn)),2);
    % ... equal to: supremum(abs(interval(VerrorDyn))); ... but faster
    
% POLY --------------------------------------------------------------------
elseif strcmp(options.alg,'poly') && options.tensorOrder == 3

    % set handles to correct files
    obj = setHessian(obj,'standard');
    obj = setThirdOrderTensor(obj,'int');
    
    % obtain intervals and combined interval z
    dz = [IH_x; IH_u];

    % compute zonotope of state and input
    Rred_diff = reduce(zonotope(Rdiff),'adaptive',sqrt(options.redFactor));
    Z_diff = cartProd(Rred_diff,options.U);
    
    % second-order error
    % use symmetry of H to simplify computation
    error_secondOrder_dyn = quadMap(Zdelta,Z_diff,H) + 0.5*quadMap(Z_diff,H);

    try
        % evaluate the third-order tensor
        if ~options.thirdOrderTensorempty
            [T,ind] = obj.thirdOrderTensor(totalInt_x, totalInt_u);
        else
            ind = {};
        end
    catch ME
        if strcmp(ME.identifier,'interval:setoutofdomain')
            ME = MException('reach:setoutofdomain',...
                'Interval Arithmetic: Set out of domain.');
            throw(ME);
        else
            rethrow(ME);
        end
    end

    % calculate the Lagrange remainder term
    if ~isempty(cell2mat(ind))
        nrind = length(ind);
        error_thirdOrder_old = interval(zeros(nrind,1),zeros(nrind,1));
        for i=1:nrind
            for j=1:length(ind{i})
                error_thirdOrder_old(i,1) = error_thirdOrder_old(i,1) + ...
                    (dz.' * T{i,ind{i}(j)} * dz) * dz(ind{i}(j));
            end
        end
        error_thirdOrder_dyn = zonotope(1/6 * error_thirdOrder_old);
    else
        error_thirdOrder_dyn = 0;
        options.thirdOrderTensorempty = true;
    end

    % combine results
    VerrorDyn = error_secondOrder_dyn + error_thirdOrder_dyn;
    VerrorDyn = reduce(VerrorDyn,'adaptive',sqrt(options.redFactor));

    temp = VerrorDyn + zonotope(VerrorStat);
    err = abs(center(temp)) + sum(abs(generators(temp)),2);
    % same as: err = supremum(abs(interval(VerrorDyn) + interval(VerrorStat)))
    
% OTHER COMBINATIONS OF ALG/TENSORORDER NOT IMPLEMENTED YET
else
    throw(CORAerror('CORA:notSupported','Specified tensor order not supported.'));
end

end


% Auxiliary functions -----------------------------------------------------

function options = aux_checkIfHessianConst(obj,options,H,totalInt_x,totalInt_u)
% check if hessian is constant --- only once executed! (very first step)
% returns: isHessianConst, hessianCheck, hessianConst

% assume true, check until false
options.isHessianConst = true;

% evaluate Hessians with larger set, check if result equal
scalingFactor = 1.1;
H_test = obj.hessian(enlarge(totalInt_x,scalingFactor),...
                     enlarge(totalInt_u,scalingFactor));
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
