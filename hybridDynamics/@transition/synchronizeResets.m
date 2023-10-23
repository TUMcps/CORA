function syncReset = synchronizeResets(transitionSet,n,m,idStates)
% synchronizeResets - merge reset functions with the same synchronization
%    labels; since the variables of the different subcomponents are
%    considered disjoint in this implementation (i.e., the left-hand side
%    of the reset function of one component is not influenced by any other
%    component), we synchronize the resets by combining all reset functions
%    into one function.
% 
%    Incoming reset functions have already been projected to the state
%    dimension of the composed automaton; thereby, all internal input
%    relations have been resolved so that the resulting reset functions
%    only depend on the full state and global inputs.
%
%    In order to obtain the correct dimensions, all reset functions need to
%    have the full state x as output dimension and the extended state [x;u]
%    composed of the full state x and the global inputs u as input
%    dimension.
%
% Syntax:
%    syncReset = synchronizeResets(transitionSet,dims,inputDims)
%
% Inputs:
%    transitionSet - array of transitions to be synchronized
%    n - full state dimension
%    m - global input dimension
%    idStates - states of components which remain as they are (identity)
%
% Outputs:
%    syncReset (struct) - synchronized reset function
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Mark Wetzlinger
% Written:       04-April-2022
% Last update:   01-July-2022
%                14-January-2023 (MW, handle states unaffected by sync)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find out which reset functions are linear and which are nonlinear
linearResets = false(size(transitionSet));
for i=1:length(transitionSet)
    if isfield(transitionSet(i).reset,'A')
        linearResets(i) = true;
    end
end
% check for dependency on inputs
withInputs = any(arrayfun(@(x) x.reset.hasInput,transitionSet,'UniformOutput',true));

% synchronize purely linear resets
if all(linearResets)

    % all reset functions need to have the same mapping
    mapDims = arrayfun(@(x) size(x.reset.A,1),transitionSet,'UniformOutput',true);
    if any(mapDims ~= mapDims(1))
        throw(CORAerror('CORA:wrongValue','first',...
            'All reset functions must map the same spaces: R^n -> R^m.'));
    end

    % initialize state/input matrix and constant offset
    A = diag(idStates);
    if withInputs
        B = zeros(n,m);
    end
    c = zeros(n,1);

    % since all resets have been projected to higher dimensions before,
    % we can just add them here
    for i=1:length(transitionSet)
        A = A + transitionSet(i).reset.A;
        if isfield(transitionSet(i).reset,'B')
            B = B + transitionSet(i).reset.B;
        end
        c = c + transitionSet(i).reset.c;
    end

    % assign matrices and vector to struct
    syncReset.A = A;
    if withInputs
        syncReset.B = B;
    end
    syncReset.c = c;
    return;
    
% convert linear resets to nonlinear resets
else
    % state dimension of all resulting nonlinear resets has to be that
    % of the extended state [x;u]
    % -> cf. evaluation in transition/reset/nonlinearReset

    for i = 1:length(linearResets)
        if linearResets(i)
            % linear resets are already projected to the full state
            % dimension -> the only remaining inputs are global inputs

            % concatenate A and B matrices for extended state vector
            if isfield(transitionSet(i).reset,'B')
                AB = [transitionSet(i).reset.A transitionSet(i).reset.B];
            else
                AB = [transitionSet(i).reset.A zeros(n,m)];
            end

            % instantiate function handle using only extended state
            % vector as input argument (note: we explicitly mention the
            % length of the extended state vector x in the anonymous
            % function to facilitate the usage of inputArgsLength in the
            % constructor call of transition below)
            if m > 0
                % at least one global input
                newReset.f = @(x,u) AB*[x(1:n);u(1:m)] + transitionSet(i).reset.c;
                newReset.hasInput = true;
            else
                % no global inputs
                newReset.f = @(x) transitionSet(i).reset.A*x(1:n) + transitionSet(i).reset.c;
                newReset.hasInput = false;
            end
            newReset.inputDim = m;
            newReset.stateDim = n;

            % instantiate nonlinear transition
            newTransition = transition(transitionSet(i).guard,newReset,...
                transitionSet(i).target,transitionSet(i).syncLabel);
            % overwrite linear transition
            transitionSet(i) = newTransition;

        end
    end
end

% synchronize nonlinear resets
    
% instantiate the output of the evaluation of the reset function
f_combined = @(x) zeros(n,1);
% add the reset functions
for j=1:length(transitionSet)
    f_combined = @(x) f_combined(x) + transitionSet(j).reset.f(x);
end
syncReset.f = @(x) f_combined(x);

% instantiate the output of the evaluation of the Jacobian
J_combined = @(x) zeros(n,n+m); 
% add the Jacobians
for i = 1:length(transitionSet)
    J_combined = @(x) J_combined(x) + transitionSet(i).reset.J(x);
end
syncReset.J = @(x) J_combined(x);

% init variable for Hessians and third-order tensor
% temp = @(x) zeros(n+m,n+m);

% Hessians
syncReset.Q = repmat({0},[n,1]);
for i = 1:n
    for j = 1:length(transitionSet)
        if ~isnumeric(transitionSet(j).reset.Q{i})
            % if Hessian is numeric, it is zero -> skip addition
            if isnumeric(syncReset.Q{i})
                % first addition to Hessian
                syncReset.Q{i} = @(x) transitionSet(j).reset.Q{i}(x);
            else
                syncReset.Q{i} = @(x) syncReset.Q{i}(x) ...
                    + transitionSet(j).reset.Q{i}(x);
            end
        end
    end
end

% third-order tensors
syncReset.T = repmat({0},[n,n+m]);
for i = 1:n
    for j = 1:n    % ...why not n+m?
        for k = 1:length(transitionSet)
            if ~isnumeric(transitionSet(k).reset.T{i,j})
                % if third-order tensor is numeric, it is zero -> skip addition
                if isnumeric(syncReset.T{i,j})
                    % first addition to third-order tensor
                    syncReset.T{i,j} = @(x) transitionSet(k).reset.T{i,j}(x);
                else
                    syncReset.T{i,j} = @(x) syncReset.T{i,j}(x) ...
                        + transitionSet(k).reset.T{i,j}(x);
                end
            end
        end
    end
end   

% properties
syncReset.stateDim = n;
syncReset.inputDim = m;
if m > 0
    syncReset.hasInput = true;
else
    syncReset.hasInput = false;
end

% ------------------------------ END OF CODE ------------------------------
