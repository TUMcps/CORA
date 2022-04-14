function obj = projectHighDim(obj,N,dim,target)
% projectHighDim - project a transition object to a higher dimensional
%                  space (required for parallel hybrid automata)
%
% Syntax:  
%    obj = projectHighDim(obj,N,dim,target,stateBinds,i)
%
% Inputs:
%    obj - transition object
%    N - dimension of the higher dimensional space
%    dim - states of the high dimensional space that correspond to the
%          states of the low dimensional transition object
%    target - target location for the new transition object

%
% Outputs:
%    obj - resulting transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, mptPolytope/projectHighDim

% Author:       Niklas Kochdumper, Johann Schoepfer
% Written:      04-January-2022
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------


    % check if the transition's reset depends on input
    if isfield(obj.reset,'hasInput') && obj.reset.hasInput
        error("Not implemented yet!");
    end

    % project guard set to the higher dimension
    guard = obj.guard;
    if ~isa(guard,'levelSet') && ~isa(guard,'mptPolytope') && ...
       ~isa(guard,'conHyperplane')
        guard = mptPolytope(guard);
    end

    obj.guard = projectHighDim(guard,N,dim);

    % project reset function to the higher dimension
    if isfield(obj.reset,'A')
        obj.reset = projectLinearReset(obj.reset,dim,N);
    else
        obj.reset = projectNonlinearReset(obj.reset,dim,N); 
    end

    % update the target location
    obj.target = target;
end


% Auxiliary Functions -----------------------------------------------------

function reset = projectLinearReset(reset,stateBind,sysDims)
% compute the projection of a linear reset function into a higher 
% dimensional space 


    % get fields from struct
    A = reset.A;
    c = reset.c;
    
    % default is identity <=> other components remain unaffected
    Aproj = eye(sysDims,sysDims);
    cProj = zeros(sysDims,1);
    
    % reset matrix A
    Aproj(stateBind,stateBind) = A;
    
    % transition vector c
    cProj(stateBind) = c;

    % construct resuling reset object
    reset.A = Aproj;
    reset.c = cProj;
    
end

function reset = projectNonlinearReset(reset,stateBinds,dims)
% compute the projection of a nonlinear reset function into a higher 
% dimensional space 

    % statebind of origin component
    stateBind = stateBinds{i};

    % project function
    reset.f = @(x) projectFunction(x,reset.f,stateBind);
    
    % project jacobian
    reset.J = @(x) projectJacobian(x,reset.J,dims,stateBind);

    % project quadratic matrix
    temp = @(x) zeros(dims);
    Q = repmat({temp},[dims,1]);
    for i = 1:length(stateBind)
       Q{stateBind(i),1} = @(x) projectQuadMat(x,reset.Q{i,1},dims, ...
                                                                stateBind); 
    end
    reset.Q = Q;
    
    % project third order tensor
    T = cell(dims,dims);
    for i = 1:length(stateBind)
        for j = 1:length(stateBind)
            if ~isempty(reset.T{i,j})
                T{stateBind(i),stateBind(j)} = @(x) projectThirdOrder(x,...
                                              reset.T{i,j},dims,stateBind);
            end
        end
    end
    reset.T = T;
end

function res = projectFunction(x,f,stateBind)
    res = x; 
    res(stateBind) = f(x(stateBind));
end

function res = projectJacobian(x,J,n,stateBind)
    res = eye(n); 
    res(stateBind,stateBind) = J(x(stateBind));
end

function res = projectQuadMat(x,Q,n,stateBind)
    res = zeros(n);
    res(stateBind,stateBind) = Q(x(stateBind));
end

function res = projectThirdOrder(x,T,n,stateBind)
    res = interval(zeros(n));
    res(stateBind,stateBind) = T(x(stateBind));
end

%------------- END OF CODE --------------