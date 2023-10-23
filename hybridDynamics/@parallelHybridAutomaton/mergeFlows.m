function sys = mergeFlows(pHA,flowList,locID)
% mergeFlows - merges the continuous dynamics of several subcomponents to 
%    obtain the continous dynamics for the overall system
%
% Syntax:
%    sys = mergeFlows(pHA,flowList,locID)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    flowList - continous dynamics object for each subcomponent
%    locID - indices of the current location
%
% Outputs:
%    sys - constructed continous dynamics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Johann Schoepfer, Niklas Kochdumper, Mark Wetzlinger
% Written:       08-June-2018  
% Last update:   09-July-2018 (NK, output instead of state for input binds)
%                22-January-2023 (MW, more general output equations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % number of components
    numComps = length(flowList);

    % check whether flow equations are linear or nonlinear
    isLinSys = false(numComps,1);
    isNonlinSys = false(numComps,1);
    for i=1:numComps
        isLinSys(i) = isa(flowList{i},'linearSys');
        isNonlinSys(i) = isa(flowList{i},'nonlinearSys');
    end

    % merge flows according to the dynamics
    if ~all(isLinSys | isNonlinSys)
        throw(CORAerror('CORA:specialError',...
            ['Only "linearSys" and "nonlinearSys" objects are currently ', ...
            'supported for parallel hybrid automata are currently supported!']));
    elseif all(isLinSys)
        sys = aux_mergeFlowsLinearSys(pHA,flowList);
    else
        % convert all linearSys to nonlinearSys
        for i=1:length(flowList)
            if isLinSys(i)
                flowList{i} = nonlinearSys(flowList{i});
            end
        end
        sys = aux_mergeFlowsNonlinearSys(pHA,flowList,locID);
    end 
end


% Auxiliary functions -----------------------------------------------------

function res = aux_mergeFlowsLinearSys(pHA,flowList)

    numComps = length(flowList);
    
    % allocate merged dynamics
    Amerged = zeros(pHA.dim,pHA.dim);
    Bmerged = zeros(pHA.dim,pHA.nrOfInputs);
    cMerged = zeros(pHA.dim,1);
    name = cell(numComps,1);
    
    % loop over all subcomponents
    for i = 1:numComps
        
        % get object properties
        flow = flowList{i};
        stateBinds = pHA.bindsStates{i};
        inputBinds = pHA.bindsInputs{i};
        
        name{i,1} = flow.name;
        A = flow.A;
        B = flow.B;
        if isempty(flow.c)
            c = zeros(size(A,1),1);
        else
            c = flow.c;
        end
        
        % constant input vector c
        cMerged(stateBinds) = cMerged(stateBinds) + c;
        
        % system matrix A
        Amerged(stateBinds,stateBinds) = A;
        
        % input matrix B
        for j = 1:size(inputBinds,1)
    
            % differentiate between global and local input
            if inputBinds(j,1) == 0
                % global input
                Bmerged(stateBinds,inputBinds(j,2)) = ...
                    Bmerged(stateBinds,inputBinds(j,2)) + B(:,j);
            
            else
                % input = output of other component, rename for clarity:
                feedComp = inputBinds(j,1);
                outputFeedComp = inputBinds(j,2);

                % flow equation of the i-th component (index 1)
                %   x1' = A1*x1 + B1*u1 + c1
                % and output equation of the feeding component (index 2)
                %   y2  = C2*x2 + D2*u2 + k2
                % with (simplified here to only one input/output)
                %   u1  = y2
                % results in the composed system
                %   x1' = A1*x1 + B1*(C2*x2 + D2*u2 + k2) + c1
                % which can be expanded to
                %   x1' = (A1*x1 + B1*C2*x2) + B1*D2*u2 + (c1 + B1*k2)
                % note: we deal with each input sequentially, so that
                % actually B1*u1 is replaced only element-wise
            
                % flow equations from the component whose output is the
                % input to the i-th component
                feedFlow = flowList{feedComp};
                % state and input binds of that component
                feedStateBinds = pHA.bindsStates{feedComp};
                feedInputBinds = pHA.bindsInputs{feedComp};
                
                % if a system has no matrices C/D/k, we assume the states
                % to be given as outputs for parallelization (just as for
                % nonlinear systems)
                if isscalar(feedFlow.C) && feedFlow.C == 1 ...
                        && isempty(feedFlow.D) && isempty(feedFlow.k)
                    C = eye(feedFlow.dim);
                else 
                    C = feedFlow.C;
                end
                
                % part with matrix C: A1*x1 + B1*C2*x2
                Amerged(stateBinds,feedStateBinds) = ...
                    Amerged(stateBinds,feedStateBinds) ...
                    + B(:,j)*C(outputFeedComp,:);
                
                % part with offset vector k: c1 + B1*k2
                if ~isempty(feedFlow.k)
                    cMerged(stateBinds) = cMerged(stateBinds) ...
                        + B(:,j)*feedFlow.k(outputFeedComp);
                end
                
                % part with feedthrough matrix D: B1*D2*u2
                if ~isempty(feedFlow.D)
                    % line in D matrix of feeding component corresponding
                    % to the bound output
                    d = feedFlow.D(outputFeedComp,:);
                    % indices of global inputs to feeding component
                    ind1 = find(feedInputBinds(:,1) == 0);
                    % check if any of the input binds of the feeding
                    % component link to other components -> not allowed!
                    ind2 = setdiff(1:size(feedInputBinds,1),ind1);
                    
                    % check D matrix regarding potential infinite loops
                    if any(d(ind2))
                        throw(CORAerror('CORA:specialError',...
                        ['It is not allowed for the feedthrough matrix D '...
                        'to point to inputs that are defined by the '...
                        'output of other subsystems, since this would '...
                        'otherwise lead to infinite loops!']));
                    end
                
                    % add to the merged B matrix from the feedthrough
                    % (binds are resolved to the global input which is fed
                    % to the feeding component)
                    if any(ind1)
                        Bmerged(stateBinds,feedInputBinds(ind1,2)) = ...
                            Bmerged(stateBinds,feedInputBinds(ind1,2)) ...
                            + B(:,j)*d(ind1);
                    end
                end
            end
        end
    end
    
    % allocate merged name (remove all default names '')
    namemerged = strjoin(name(cellfun(@(x)~isempty(x),name,'UniformOutput',true)),' x ');
    
    % construct resulting continuous dynamics object
    res = linearSys(namemerged,Amerged,Bmerged,cMerged);
end

function res = aux_mergeFlowsNonlinearSys(pHA,flowList,locID)

    % number of components in parallel hybrid automaton
    numComps = length(flowList);
    
    % construct symbolic state vector and input vector for merged flow
    x = sym('x',[pHA.dim,1]);
    u = sym('u',[pHA.nrOfInputs,1]);
    
    % initialize dynamic function
    f = sym(zeros(pHA.dim,1));
    
    % loop over all subcomponents
    for i = 1:numComps
    
        % get object properties
        flow = flowList{i};
        stateBinds = pHA.bindsStates{i};
        inputBinds = pHA.bindsInputs{i};
        
        % construct input vector for this subcomponent
        u_ = sym(zeros(size(inputBinds,1),1));
        
        % resolve inputs to either global input or output equation of
        % location of another component
        for j = 1:size(inputBinds,1)
        
            % component from which the input is coming from (0 = global)
            feedingComp = inputBinds(j,1);

            if feedingComp == 0
                % global input
                u_(j) = u(inputBinds(j,2));

            else
                % input to be taken from output equation of other component
                
                % flow equation of feeding component (correct location)
                feedingFlow = flowList{feedingComp};
                
                % state binds from feeding component
                feedingCompStateBinds = pHA.bindsStates{inputBinds(j,1)};
                % symbolic variables for states (resolved)
                xTemp = x(feedingCompStateBinds);
                % symbolic variables for inputs (unresolved)
                uTemp = sym('u',[feedingFlow.nrOfInputs,1]);

                % insert symbolic variables into output equation
                outTemp = feedingFlow.out_mFile(xTemp,uTemp);
                % output equation corresponding to input
                outTemp = outTemp(inputBinds(j,2));

                % check dependency on (local) inputs
                if logical(outTemp - subs(outTemp,uTemp,zeros(feedingFlow.nrOfInputs,1)) == 0)
                    % output equation does not depend on inputs
                    u_(j) = outTemp;

                else

                    % resolve inputs from output equation of feeding component
                    feedingCompInputBinds = pHA.bindsInputs{inputBinds(j,1)};

                    % loop over all inputs
                    for uu=1:feedingFlow.nrOfInputs
                        % check if uu-th input used in output equation
                        if logical(outTemp - subs(outTemp,uTemp(uu),0) ~= 0)
                            if feedingCompInputBinds(uu,1) == 0
                                % global input
                                outTemp = subs(outTemp,uTemp(uu),...
                                    u(feedingCompInputBinds(uu,2)));
    
                            else
                                % local input = output of another component
                                % ...currently not supported (would need to
                                % be traced back until everything resolves
                                % to states or global inputs... caution:
                                % there is a potential for infinite loops)
                                throw(CORAerror('CORA:notSupported',...
                                    ['In a parallel hybrid automaton, '...
                                     'the output equations of the flow '...
                                     'in each location can currently '...
                                     'only depend on global inputs.']));
    
                            end
                        end

                    end

                end

            end
        end
        
        % construct flow function for this subcomponent
        f(stateBinds,:) = flow.mFile(x(stateBinds),u_);
    end
    
    % name of resulting nonlinearSys object (note: we cannot use the same
    % construction as for linear systems (e.g., 'sys1 x sys2'), because the
    % name is used for the function handle -> must not contain any spaces
    name = ['location',strrep(strcat(num2str(locID')),' ','')];
    % function handle for flow equation
    funHan = matlabFunction(f,'Vars',{x,u});
    
    % instantiate nonlinear system
    res = nonlinearSys(name,funHan,pHA.dim,pHA.nrOfInputs);
end

% ------------------------------ END OF CODE ------------------------------
