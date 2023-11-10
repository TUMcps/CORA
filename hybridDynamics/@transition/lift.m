function trans = lift(trans,N,dims,target,id)
% lift - project a transition object to a higher dimensional space
%    (required for parallel hybrid automata); the reset function of the
%    transition does not depend on inputs (e.g., states of other components)
%
% Syntax:
%    trans = lift(trans,N,dims,target,id)
%
% Inputs:
%    trans - transition object
%    N - dimension of the higher dimensional space
%    dims - states of the high dimensional space that correspond to the
%           states of the low dimensional transition object
%    target - target location for the new transition object
%    id - true/false whether identity reset function should be used for all
%         other states
%
% Outputs:
%    trans - resulting transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, polytope/lift,
%           projectInputDependentTrans

% Authors:       Niklas Kochdumper, Johann Schoepfer
% Written:       04-January-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check if the transition's reset depends on input
    if trans.reset.hasInput
        throw(CORAerror('CORA:specialError',...
            ['The projection of input-dependent transitions is only '...
            'implemented within parallel hybrid automata.']));
        % The projection of input-dependent transitions is implemented for
        % use within pHAs (see projectInputDependentTrans.m), but not for
        % use on a transition alone, because the current implementation
        % is specialized on handling inputs coming from other components
        % from within the pHA).
    end

    
    % convert to polytope unless guard set is fullspace, a levelSet,
    % a polytope, or a conHyperplane
    guard = trans.guard;
    if ~isa(guard,'fullspace') && ~isa(guard,'levelSet') ...
            && ~isa(guard,'polytope') && ~isa(guard,'conHyperplane')
        guard = polytope(guard);
    end
    
    % project guard set to the higher dimension
    trans.guard = lift(guard,N,dims);

    % project reset function to the higher dimension
    if isfield(trans.reset,'A')
        trans.reset = aux_liftLinearReset(trans.reset,dims,N,id);
    else
        trans.reset = aux_liftNonlinearReset(trans.reset,dims,N,id); 
    end

    % update the target location
    trans.target = target;
end


% Auxiliary functions -----------------------------------------------------

function resetStruct = aux_liftLinearReset(resetStruct,stateBind,dims,id)
% lits a linear reset function into a higher-dimensional space

    % note: reset function needs to map from R^n -> R^n
    [n_,n] = size(resetStruct.A);
    if n_ ~= n
        throw(CORAerror('CORA:notSupported',...
            ['Projection of reset functions to higher-dimensional spaces '...
            'only supported for R^n -> R^n.']));
    end
    
    % init A matrix by identity or zeros
    if id
        Aproj = eye(dims,dims);
    else
        Aproj = zeros(dims,dims);
    end
    cProj = zeros(dims,1);
    
    % insert reset matrix A and vector c
    Aproj(stateBind,stateBind) = resetStruct.A;
    cProj(stateBind) = resetStruct.c;

    % construct resulting reset object
    resetStruct.A = Aproj;
    resetStruct.c = cProj;

    % update properties
    resetStruct.stateDim = dims;
    resetStruct.inputDim = 0;
    resetStruct.hasInput = false;
end

function resetStruct = aux_liftNonlinearReset(resetStruct,stateBind,dims,id)
% lifts a nonlinear reset function into a higher-dimensional space
% note: currently no option to set other states to identity?

    % projection only allowed for R^n -> R^n (and no inputs)
    [in,out] = inputArgsLength(resetStruct.f);
    if in(2) > 0
        throw(CORAerror('CORA:notSupported',...
            ['Projection of reset functions to higher-dimensional spaces '...
            'only for state-dependent resets']));
    elseif in(1) ~= out
        throw(CORAerror('CORA:notSupported',...
            ['Projection of reset functions to higher-dimensional spaces '...
            'only supported for R^n -> R^n.']));
    end

    % project function
    resetStruct.f = @(x) aux_projectFunction(x,resetStruct.f,dims,stateBind);
    
    % project jacobian
    resetStruct.J = @(x) aux_projectJacobian(x,resetStruct.J,dims,stateBind);

    % project hessian matrix
    temp = @(x) zeros(dims);
    Q = repmat({temp},[dims,1]);
    for i = 1:length(stateBind)
        Q{stateBind(i),1} = ...
            @(x) aux_projectQuadMat(x,resetStruct.Q{i,1},dims,stateBind); 
    end
    resetStruct.Q = Q;
    
    % project third-order tensor
    T = cell(dims,dims);
    for i = 1:length(stateBind)
        for j = 1:length(stateBind)
            if ~isempty(resetStruct.T{i,j})
                T{stateBind(i),stateBind(j)} = ...
                    @(x) aux_projectThirdOrder(x,resetStruct.T{i,j},dims,stateBind);
            end
        end
    end
    resetStruct.T = T;

    % update properties
    resetStruct.stateDim = dims;
    resetStruct.inputDim = 0;
    resetStruct.hasInput = false;
end

function res = aux_projectFunction(x,f,n,stateBind)
% project reset function to higher dimension
    res = x(1:n); 
    res(stateBind) = f(x(stateBind));
end

function res = aux_projectJacobian(x,J,n,stateBind)
% project Jacobian of reset function to higher dimension
    res = eye(n); 
    res(stateBind,stateBind) = J(x(stateBind));
end

function res = aux_projectQuadMat(x,Q,n,stateBind)
% project Hessian of reset function to higher dimension
    res = zeros(n);
    try
        res(stateBind,stateBind) = Q(x(stateBind));
    end
end

function res = aux_projectThirdOrder(x,T,n,stateBind)
% project third-order tensor reset function to higher dimension
    res = interval(zeros(n));
    try
        res(stateBind,stateBind) = T(x(stateBind));
    end
end

% ------------------------------ END OF CODE ------------------------------
