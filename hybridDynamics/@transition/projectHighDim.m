function trans = projectHighDim(trans,N,dims,target,id)
% projectHighDim - project a transition object to a higher dimensional
%    space (required for parallel hybrid automata); the reset function of
%    the transition does not depend on inputs (e.g., states of other
%    components)
%
% Syntax:  
%    trans = projectHighDim(trans,N,dims,target,id)
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
% See also: parallelHybridAutomaton, mptPolytope/projectHighDim,
%           projectInputDependentTrans

% Author:       Niklas Kochdumper, Johann Schoepfer
% Written:      04-January-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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

    % project guard set to the higher dimension
    guard = trans.guard;
    if ~(isnumeric(guard) && isempty(guard)) && ~isa(guard,'levelSet') ...
            && ~isa(guard,'mptPolytope') && ~isa(guard,'conHyperplane')
        % convert to mptPolytope unless guard set is [], a levelSet,
        % an mptPolytope, or a conHyperplane
        guard = mptPolytope(guard);
    end
    
    % only project unless guard set is []
    if ~(isnumeric(guard) && isempty(guard))
        trans.guard = projectHighDim(guard,N,dims);
    end

    % project reset function to the higher dimension
    if isfield(trans.reset,'A')
        trans.reset = projectLinearReset(trans.reset,dims,N,id);
    else
        trans.reset = projectNonlinearReset(trans.reset,dims,N,id); 
    end

    % update the target location
    trans.target = target;
end


% Auxiliary Functions -----------------------------------------------------

function resetStruct = projectLinearReset(resetStruct,stateBind,dims,id)
% projects a linear reset function into a higher-dimensional space 
    
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

function resetStruct = projectNonlinearReset(resetStruct,stateBind,dims,id)
% compute the projection of a nonlinear reset function into a higher 
% dimensional space
% note: currently no option to set other states to identity?

    % project function
    resetStruct.f = @(x) projectFunction(x,resetStruct.f,stateBind);
    
    % project jacobian
    resetStruct.J = @(x) projectJacobian(x,resetStruct.J,dims,stateBind);

    % project hessian matrix
    temp = @(x) zeros(dims);
    Q = repmat({temp},[dims,1]);
    for i = 1:length(stateBind)
        Q{stateBind(i),1} = ...
            @(x) projectQuadMat(x,resetStruct.Q{i,1},dims,stateBind); 
    end
    resetStruct.Q = Q;
    
    % project third-order tensor
    T = cell(dims,dims);
    for i = 1:length(stateBind)
        for j = 1:length(stateBind)
            if ~isempty(resetStruct.T{i,j})
                T{stateBind(i),stateBind(j)} = ...
                    @(x) projectThirdOrder(x,resetStruct.T{i,j},dims,stateBind);
            end
        end
    end
    resetStruct.T = T;

    % update properties
    resetStruct.stateDim = dims;
    resetStruct.inputDim = 0;
    resetStruct.hasInput = false;
end

function res = projectFunction(x,f,stateBind)
% project reset function to higher dimension
    res = x; 
    res(stateBind) = f(x(stateBind));
end

function res = projectJacobian(x,J,n,stateBind)
% project Jacobian of reset function to higher dimension
    res = eye(n); 
    res(stateBind,stateBind) = J(x(stateBind));
end

function res = projectQuadMat(x,Q,n,stateBind)
% project Hessian of reset function to higher dimension
    res = zeros(n);
    res(stateBind,stateBind) = Q(x(stateBind));
end

function res = projectThirdOrder(x,T,n,stateBind)
% project third-order tensor reset function to higher dimension
    res = interval(zeros(n));
    res(stateBind,stateBind) = T(x(stateBind));
end

%------------- END OF CODE --------------