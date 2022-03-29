function trans = projectInputDependentTrans(pha, trans,N,compIndex,loc,targ)
% projectInputDependentTrans - Constructs a transition of adequat dimension
%                              for a PHA with a transition which depends on
%                              states and inputs
%
% Syntax:
%    res = projectInputDependentTrans(pha, transition)
%
% Inputs:
%    pha - parallel hybrid automaton object
%    transition - the transition in question
%    N - dimensions of the higher dimensional space
%    compIndex - index of the component belonging to the transition
%    loc - location vector of the pha
%    targ - target of the generated transition
%    
%
% Outputs:
%    transition - constructed transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Maximilian Perschl
% Written:      17-February-2022
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % project guard set to the higher dimension
    guard = trans.guard;
    if ~isa(guard,'levelSet') && ~isa(guard,'mptPolytope') && ...
            ~isa(guard,'conHyperplane')
        guard = mptPolytope(guard);
    end

    guard = projectHighDim(guard,N,pha.bindsStates{compIndex});

    % project reset function to the higher dimension
    if isfield(trans.reset,'A')
        reset = projectLinearReset(pha,trans.reset,compIndex,loc);
    else
        reset = projectNonlinearReset(pha,trans.reset,compIndex);
    end
    
    % update the target location
   target = targ;
   
   % create transition object
   trans = transition(guard,reset,target);
   
end

% Auxiliary Functions -----------------------------------------------------
function reset = projectLinearReset(pha,reset,compIndex,loc)
% compute new reset matrices A,B,c according to the state/input dimensions
% of the new automata, with the mapping of outputs of other components to
% the inputs of the automaton which owns the relevant transition taken into
% account

    % get fields from struct
    A = reset.A;
    B = reset.B;
    c = reset.c;
    
    stateBinds = pha.bindsStates;
    inputBinds = pha.bindsInputs;
    
    sysDim = pha.numStates;
    sysInput = pha.numInputs;
    
    % Instantiate new matrices
    % default is identity <=> other components remain unaffected
    Aproj = eye(sysDim,sysDim);
    Bproj = zeros(sysDim,sysInput);
    cProj = zeros(sysDim,1);
    
    % define reset matrix A
    % old matrix A is taken as is, just at the correct indices concerning
    % stateBinds
    Aproj(stateBinds{compIndex},stateBinds{compIndex}) = A;
    % sort inputs by global and received from other components
    bind = inputBinds{compIndex};
    for i = 1:size(bind,1)
        if bind(i,1) == 0
            % global input --> old B matrix in new B matrix
            globalIndex = bind(i,2);
            Bproj(stateBinds{compIndex},globalIndex) = B(:,i);
        else
            % input received from subcomponent --> set entries in Aproj and
            % Bproj equal to subcomponent's flow's C and D multiplied by
            % transition's B
            outsideComponentIndex = bind(i,1);
            outsideOutputIndex = bind(i,2);
            % get dynamics of other subcomponent
            outsideDynamics = pha.components{outsideComponentIndex}.location{loc(outsideOutputIndex)}.contDynamics;
            % other subcomponent flow linear or nonlinear?
            if isa(outsideDynamics,'linearSys')
                % reset is Ax + Bu + c
                % u is C'*x + D'*u' + k' with C',D',k',u' being 
                % output dynamics and inputs to outside component
                % => Bu = BC'x + BD'u' + Bk'
                
                
                % ------Note: this check might be irrelevant because-------------------
                %             it is already checked in mergeFlows.m
                
                % if D'u' contains non-global inputs we give an error, since this
                % could lead to infinite loops 
                % (either D'(_,i) needs to be zero or u'(i) needs to be
                % global)
                for j = 1:size(outsideDynamics.D,2)
                    if ~(outsideDynamics.D(outsideOutputIndex,j) == 0 || inputBinds{outsideComponentIndex}(j,1) == 0)
                        error(['It is not allowed for the throughput matrix D '...
                            'to point to inputs that are defined by the '...
                            'output of other subsystems, since it would '...
                            'otherwise be able to construct infinite loops!']);                 
                    end
                end
                outsideSystem_C = outsideDynamics.C; 
                outsideSystem_D = outsideDynamics.D;
                outsideSystem_k = outsideDynamics.k;
                
                % if a system has no matrices C/D/k, we assume the states are
                % given as outputs for parallelization (just as for nonlinear
                % systems)
                if outsideSystem_C == 1 && isempty(outsideSystem_D) && isempty(outsideSystem_k)
                    outsideSystem_C = eye(outsideDynamics.dim);
                end
                % compute adjusted matrices (only one column needed)
                BC_i = B(:,i).*outsideSystem_C(outsideOutputIndex,:);
                BD_i = 0;
                if ~isempty(outsideSystem_D)
                	BD_i = B(:,i).*outsideSystem_D(outsideOutputIndex,:);
                end
                Bk_i = 0;
                if ~isempty(outsideSystem_k)
                    Bk_i = B(:,i).*outsideSystem_k(outsideOutputIndex);
                end
                % adjust global reset matrices
                
                Aproj(stateBinds{compIndex},stateBinds{outsideComponentIndex}) = ...
                    Aproj(stateBinds{compIndex},stateBinds{outsideComponentIndex}) + BC_i;
                
                
                % get indices of global inputs outside component receives 
                globalInputIndices = [];
                for j = 1:size(inputBinds{outsideComponentIndex},1)
                    if inputBinds{outsideComponentIndex}(j,1) == 0
                        globalInputIndices = vertcat(globalInputIndices,inputBinds{outsideComponentIndex}(j,2));
                    end
                end
                if ~isempty(globalInputIndices)
                    Bproj(stateBinds{compIndex},inputBinds{outsideComponentIndex}) = ...
                        Bproj(stateBinds{compIndex},globalInputIndices) + BD_i;
                end
                cProj(stateBinds{compIndex}) = cProj(stateBinds{compIndex}) + Bk_i;
            else
                error("Mixing different system classes is not allowed for PHAs!");
            end
        end
    end
    
    
    % transition vector c
    cProj(stateBinds{compIndex}) = c;
    
    % update reset struct
    reset.A = Aproj;
    reset.B = Bproj;
    reset.c = cProj;
    reset.hasInput = 1;
    reset.inputDim = pha.numInputs;
end


function resetResult = projectNonlinearReset(pha,reset,compIndex)
% project nonlinear reset function dependent on inputs
% change reset function along with derivatives up to third order 

% compute state and input variables
x = sym('x',[pha.numStates,1]);
u = sym('u',[pha.numInputs,1]);

% initialize dynamic function
f = sym(ones(pha.numStates,1));

% compute symbolic vector resolving input mapping
% vector models the altered state x' = [x;u] of the subcomponent
x_reset = x(pha.bindsStates{compIndex});

% add inputs/outside states
% Note: Nonlinear-systems do not have outputs right now, so we take the
% i-th state instead
for i = 1:size(pha.bindsInputs{compIndex},1)
    binds = pha.bindsInputs{compIndex}(i,:);
    if binds(1) == 0 
        x_reset = [x_reset;u(binds(2))];
    else
        x_reset = [x_reset;x(pha.bindsStates{binds(1)}(binds(2)))];
    end
end


% get symbolic function of original reset function
% with inputs mapping resolved
fOld = reset.f(x_reset);

% inject original reset in projected reset
f(pha.bindsStates{compIndex}) = fOld;

% convert symbolic function to function handle,
f_han = matlabFunction(f,'Vars',{x,u});

% compute jacobian, hessian and third order tensor
% with augmented state x' = [x;u];

x_aug = [x;u];

J = jacobian(f,x_aug);

J_han = matlabFunction(J,'Vars',{x_aug});

Q = cell(pha.numStates,1);

for i = 1:pha.numStates
    Q{i,1} = matlabFunction(hessian(f(i),x),'Vars',{x_aug});
end

T = cell(size(J,1));

for i = 1:size(J,1)
    for j = 1:size(J,2)
        temp = hessian(J(i,j),x_aug);
        if any(any(temp~=0))
            T{i,j} = matlabFunction(temp,'Vars',{x_aug});
        end
    end
end

% create resulting reset struct
resetResult.f = f_han;

resetResult.J = J_han;

resetResult.Q = Q;

resetResult.T = T;

resetResult.hasInput = 1;

resetResult.inputDim = pha.numInputs;

end
%------------- END OF CODE --------------