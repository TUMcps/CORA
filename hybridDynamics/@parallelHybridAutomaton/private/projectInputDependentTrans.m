function trans = projectInputDependentTrans(pHA,trans,N,compIndex,locID,targ,id)
% projectInputDependentTrans - Constructs a transition of adequate
%    dimension for a parallel hybrid automaton with a transition which
%    depends on inputs (e.g., states of other components)
%
% Syntax:
%    res = projectInputDependentTrans(pHA,trans,N,compIndex,locID,targ,id)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    trans - the transition in question
%    N - dimensions of the higher dimensional space
%    compIndex - index of the component belonging to the transition
%    locID - IDs of the currently active locations
%    targ - target of the generated transition
%    id - true/false whether identity reset function should be used for all
%         other states
%
% Outputs:
%    transition - constructed transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       17-February-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % convert to polytope unless guard set is fullspace, a levelSet,
    % a polytope, or a conHyperplane
    guard = trans.guard;
    if ~isa(guard,'fullspace') && ~isa(guard,'levelSet') ...
            && ~isa(guard,'polytope') && ~isa(guard,'conHyperplane')
        guard = polytope(guard);
    end

    % project guard set to the higher dimension
    guard = lift_(guard,N,pHA.bindsStates{compIndex});

    % project reset function to the higher dimension
    if isfield(trans.reset,'A')
        reset = aux_projectLinearReset(pHA,trans.reset,compIndex,locID,id);
    else
        reset = aux_projectNonlinearReset(pHA,trans.reset,compIndex,locID,id);
    end
    
    % update the target location
    target = targ;
   
    % create transition object
    trans = transition(guard,reset,target);
   
end


% Auxiliary functions -----------------------------------------------------

function resetStruct = aux_projectLinearReset(pHA,resetStruct,compIndex,locID,id)
% compute new reset matrices A,B,c according to the state/input dimensions
% of the new automaton, where the mapping of inputs (which are states of
% other components, defined in pHA.bindsInputs) is resolved; the resulting
% reset function depends only on the full state and the global inputs

    % get fields from struct
    A = resetStruct.A;
    B = resetStruct.B;
    c = resetStruct.c;
    
    % read state and input binds
    stateBinds = pHA.bindsStates;
    inputBinds = pHA.bindsInputs;
    
    % read number of states and inputs of the entire automaton
    sysDim = pHA.dim;
    sysInput = pHA.nrOfInputs;
    
    % Instantiate matrices in full state dimension
    if id
        Aproj = eye(sysDim,sysDim);
    else
        Aproj = zeros(sysDim,sysDim);
    end
    Bproj = zeros(sysDim,sysInput);
    cProj = zeros(sysDim,1);
    
    % define reset matrix A
    % old matrix A is taken as is, just at the correct indices concerning
    % stateBinds
    Aproj(stateBinds{compIndex},stateBinds{compIndex}) = A;
    % sort inputs by global and received from other components
    bind = inputBinds{compIndex};

    % loop over all inputs to current component
    for i = 1:size(bind,1)
        if bind(i,1) == 0
            % input i is a global input --> insert that part of the old
            % input matrix into the new input matrix
            Bproj(stateBinds{compIndex},bind(i,2)) = B(:,i);

        else
            % input i is the state of another component --> set entries in
            % Aproj and Bproj equal to that component's matrices C and D
            % multiplied by the transition's B

            % index of component where input i comes from
            otherCompIdx = bind(i,1);
            % row in output equation of other component
            otherOutputIdx = bind(i,2);

            % get flow of other component
            otherFlow = pHA.components(otherCompIdx).location(locID(otherCompIdx)).contDynamics;

            % is other component's flow linear or nonlinear?
            if isa(otherFlow,'linearSys')
                % our reset is
                %   Ax + Bu + c
                % where the input u is
                %   C'*x + D'*u' + k'
                % where C',D',k',u' are variables of the other component;
                % the resulting value for Bu in our reset then becomes:
                %   Bu = BC'x + BD'u' + Bk'
                
                % if D'u' contains non-global inputs we give an error,
                % since this could lead to infinite loops 
                % (either D'(_,i) needs to be zero or u'(i) needs to be
                % global)
                for j = 1:size(otherFlow.D,2)
                    if ~(otherFlow.D(otherOutputIdx,j) == 0 ...
                            || inputBinds{otherCompIdx}(j,1) == 0)
                        throw(CORAerror('CORA:notSupported',...
                            ['It is not allowed for the feedthrough matrix D '...
                            'to point to inputs that are defined by the '...
                            'output of other subsystems, since it would '...
                            'otherwise be able to construct infinite loops!']));                 
                    end
                end

                % read out output equation of other component
                otherFlow_C = otherFlow.C; 
                otherFlow_D = otherFlow.D;
                otherFlow_k = otherFlow.k;
                
                % if a system has no matrices C/D/k, we assume the states
                % are given as outputs for parallelization (just as for
                % nonlinear systems)
                if isscalar(otherFlow_C) && otherFlow_C == 1 ...
                        && isempty(otherFlow_D) && isempty(otherFlow_k)
                    otherFlow_C = eye(otherFlow.dim);
                end

                % compute adjusted matrices (only requires the column
                % corresponding to the state our reset function receives as
                % an input)
                BC_i = B(:,i) .* otherFlow_C(otherOutputIdx,:);
                BD_i = 0;
                if ~isempty(otherFlow_D)
                	BD_i = B(:,i) .* otherFlow_D(otherOutputIdx,:);
                end
                Bk_i = 0;
                if ~isempty(otherFlow_k)
                    Bk_i = B(:,i) .* otherFlow_k(otherOutputIdx);
                end

                % adjust full reset matrices A, B, and vector c
                Aproj(stateBinds{compIndex},stateBinds{otherCompIdx}) = ...
                    Aproj(stateBinds{compIndex},stateBinds{otherCompIdx}) + BC_i;
                
                % indices of the global inputs the other component receives 
                globalInputIndices = [];
                for j = 1:size(inputBinds{otherCompIdx},1)
                    if inputBinds{otherCompIdx}(j,1) == 0
                        globalInputIndices = vertcat(globalInputIndices,...
                            inputBinds{otherCompIdx}(j,2));
                    end
                end
                if ~isempty(globalInputIndices)
                    Bproj(stateBinds{compIndex},globalInputIndices) = ...
                        Bproj(stateBinds{compIndex},globalInputIndices) ...
                        + BD_i(:,globalInputIndices);
                end
                cProj(stateBinds{compIndex}) = cProj(stateBinds{compIndex}) + Bk_i;

            else
                throw(CORAerror('CORA:notSupported',...
                    ['Mixing different system classes is not supported for '...
                    'reachability analysis of parallel hybrid automata.']));
            end

        end
    end
    
    % transition vector c
    cProj(stateBinds{compIndex}) = c;
    
    % assign matrices and vector
    resetStruct.A = Aproj;
    resetStruct.B = Bproj;
    resetStruct.c = cProj;

    % update properties
    resetStruct.stateDim = sysDim;    
    resetStruct.inputDim = sysInput;
    if sysInput > 0
        resetStruct.hasInput = true;
    else
        resetStruct.hasInput = false;
    end

end

function resetResult = aux_projectNonlinearReset(pHA,resetStruct,compIdx,locID,id)
% project nonlinear reset function to the full state dimension, where the
% reset function depends on inputs (either states of other components or
% global inputs); additionally, we update all derivatives up to third order

    % enlist state and input variables
    x = sym('x',[pHA.dim,1]);
    u = sym('u',[pHA.nrOfInputs,1]);
    
    % initialize reset function
    if id
        f = x;
    else
        f = sym(zeros(pHA.dim,1));
    end

    % read out state binds
    statebinds = pHA.bindsStates{compIdx};
    
    % initialize symbolic vector containing of states of given component
    x_reset = x(statebinds);
    
    % resolve input binds: in addition to its own states, the reset
    % function of the component takes the inputs (either outputs from other
    % components or global inputs) as input arguments
    for i = 1:size(pHA.bindsInputs{compIdx},1)
        % read out input bind:
        % binds(1) is component from which the output is fed as an input
        % binds(2) is the number of the fed output from that component
        binds = pHA.bindsInputs{compIdx}(i,:);
        
        if binds(1) == 0
            % global input
            x_reset = [x_reset;u(binds(2))];
        else
            % output from another component: insert state vector into
            % corresponding output equation
            othersys = pHA.components(binds(1)).location(locID(binds(1))).contDynamics;
            if isa(othersys,'linearSys')
                y = othersys.C * x(pHA.bindsStates{binds(1)}) ...
                    + othersys.D * u;
            elseif isa(othersys,'nonlinearSys')
                y = othersys.out_mFile(x(pHA.bindsStates{binds(1)}),u);
            end

            % concatenate
            % x_reset = [x_reset;x(pHA.bindsStates{binds(1)}(binds(2)))];
            x_reset = [x_reset;y(binds(2))];
        end
    end
    
    % symbolic evaluation of the component's reset function (input mappings
    % are resolved and now corresponds to state variables)
    fOld = resetStruct.f(x_reset);
    
    % insert original reset in projected reset
    f(statebinds) = fOld;
    
    % augment the state: full state + global inputs
    x_aug = [x;u];

    % convert symbolic function to function handle
    f_han = matlabFunction(f,'Vars',{x_aug});
    
    % compute Jacobian
    J = jacobian(f,x_aug);
    J_han = matlabFunction(J,'Vars',{x_aug});
    
    % compute Hessian
    Q = cell(pHA.dim,1);
    % skip all dimensions that do not occur in state binds (= all-zero)
    for i=1:length(statebinds)
        temp = hessian(f(statebinds(i)),x_aug);
        Q{statebinds(i),1} = matlabFunction(temp,'Vars',{x_aug});
    end
    
    % compute third-order tensor
    T = cell(pHA.dim,pHA.dim+pHA.nrOfInputs);
    % skip all rows in Jacobian that do not occur in state binds
    for i = 1:length(statebinds)
        for j = 1:pHA.dim+pHA.nrOfInputs
            temp = hessian(J(statebinds(i),j),x_aug);
            if any(any(temp~=0))
                T{statebinds(i),j} = matlabFunction(temp,'Vars',{x_aug});
            end
        end
    end
    
    % instantiate reset struct
    resetResult.f = f_han;
    resetResult.J = J_han;
    resetResult.Q = Q;
    resetResult.T = T;

    % properties
    resetResult.stateDim = pHA.dim;
    resetResult.inputDim = pHA.nrOfInputs;
    if pHA.nrOfInputs > 0
        resetResult.hasInput = true;
    else
        resetResult.hasInput = false;
    end

end

% ------------------------------ END OF CODE ------------------------------
