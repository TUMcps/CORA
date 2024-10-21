function sys = mergeFlows(pHA,locID)
% mergeFlows - merges the continuous dynamics of several subcomponents to 
%    obtain the continous dynamics for the overall system
%
% Syntax:
%    sys = mergeFlows(pHA,locID)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
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
    numComps = length(pHA.components);

    % read out flow of active location of i-th subcomponent
    flowList = cell(1,numComps);
    for i=1:numComps
        flowList{i} = pHA.components(i).location(locID(i)).contDynamics;
    end

    % check whether flow equations are linear or nonlinear
    isLinsys = false(numComps,1);
    isNonlinsys = false(numComps,1);
    for i=1:numComps
        isLinsys(i) = isa(flowList{i},'linearSys');
        isNonlinsys(i) = isa(flowList{i},'nonlinearSys');
    end

    % merge flows according to the dynamics
    if ~all(isLinsys | isNonlinsys)
        throw(CORAerror('CORA:specialError',...
            ['Only "linearSys" and "nonlinearSys" objects are currently ', ...
            'supported for parallel hybrid automata are currently supported!']));
    elseif all(isLinsys)
        sys = aux_mergeFlowsLinearSys(pHA,flowList);
    else
        % convert all linearSys to nonlinearSys
        for i=1:length(flowList)
            if isLinsys(i)
                flowList{i} = nonlinearSys(flowList{i});
            end
        end
        sys = aux_mergeFlowsNonlinearSys(pHA,flowList,locID);
    end 
end


% Auxiliary functions -----------------------------------------------------

function sys_merged = aux_mergeFlowsLinearSys(pHA,flowList)
% merge the flows of all linear systems from the current location ID to a
% joint flow; this becomes difficult because inputs of some components are
% outputs of others
% note: disturbances and noises are globally unique per component, that is,
% disturbance i to component j is only used there and in no other component

    numComps = length(flowList);
    
    % allocate merged dynamics
    name = cell(numComps,1);
    Amerged = zeros(pHA.nrOfStates,pHA.nrOfStates);
    Bmerged = zeros(pHA.nrOfStates,pHA.nrOfInputs);
    cMerged = zeros(pHA.nrOfStates,1);
    Emerged = zeros(pHA.nrOfStates,pHA.nrOfDisturbances);
    Fmerged = zeros(pHA.nrOfStates,pHA.nrOfNoises);

    % save index for disturbances
    idxDist = 1;
    % pre-compute noise binds
    cumulativeDists = [0; ...
        cumsum(arrayfun(@(x) x.location(1).contDynamics.nrOfDisturbances, pHA.components))];
    bindsNoises = cell(numComps,1);
    for i = 1:numComps
        bindsNoises{i} = (cumulativeDists(i)+1):cumulativeDists(i+1);
    end
    
    % loop over all subcomponents
    for i = 1:numComps
        
        % get object properties
        flow = flowList{i};
        stateBinds = pHA.bindsStates{i};
        inputBinds = pHA.bindsInputs{i};
        
        name{i,1} = flow.name;
        A = flow.A;
        B = flow.B;
        c = flow.c;
        E = flow.E;
        
        % constant input vector c
        cMerged(stateBinds) = cMerged(stateBinds) + c;

        % system matrix A
        Amerged(stateBinds,stateBinds) = A;

        % disturbance matrix E
        distComp = size(E,2);
        Emerged(stateBinds,idxDist:idxDist+distComp-1) = E;
        idxDist = idxDist + distComp;
        
        % input binds: system matrix A (via C), input matrix B, offset c
        for j = 1:size(inputBinds,1)
    
            % distinguish between global and local input
            if inputBinds(j,1) == 0
                % global input
                Bmerged(stateBinds,inputBinds(j,2)) = ...
                    Bmerged(stateBinds,inputBinds(j,2)) + B(:,j);
            
            else
                % input = output of other component, rename for clarity:
                feedComp = inputBinds(j,1);
                outputFeedComp = inputBinds(j,2);

                % flow equation of the i-th component (index 1)
                %   x1'(i) = A1(i,:)*x1 + B1(i,:)*u1 + c1(i) + E1(i,:)*w1
                % and output equation of the feeding component (index 2)
                %   y2(k)  = C2(k,:)*x2 + D2(k,:)*u2 + k2(k) + F2(k,:)*v2
                % with -- looking here to only one input/output pair --
                %   u1(j) = y2(k)
                % results in the composed system
                %   x1'(i) = A1(i,:)*x1
                %            + B1(i,:) * (C2(k,:)*x2 + D2(k,:)*u2 + k2 + F2(k,:)*v2)
                %            + c1(i) + E1(i,:)*w1
                % which can be expanded to
                %   x1'(i) = (A1(i,:)*x1 + B1(i,:)*C2(k,:)*x2)
                %            + B1(i,:)*D2(k,:)*u2
                %            + (c1(i) + B1(i,:)*k2(k))
                %            + E1(i,:)*w1
                %            + B1(i,:)*F2(k,:)*v2
                % note: we deal with each input sequentially, so that
                % actually B1*u1 is replaced only element-wise (iterator: j)
            
                % flow equations from the component whose output is the
                % input to the i-th component
                feedFlow = flowList{feedComp};
                % state, input, and noise binds of that component
                feedStateBinds = pHA.bindsStates{feedComp};
                feedInputBinds = pHA.bindsInputs{feedComp};
                feedNoiseBinds = bindsNoises{feedComp};
                
                % if a system has no matrices D/k/F, we assume the states
                % to be given as outputs for parallelization (just as for
                % nonlinear systems)
                if isscalar(feedFlow.C) && feedFlow.C == 1 ...
                        && ~any(any(feedFlow.D)) && ~any(feedFlow.k) ...
                        && ~any(any(feedFlow.F))
                    C = eye(feedFlow.nrOfStates);
                else 
                    C = feedFlow.C;
                end
                
                % part with matrix C:
                %    A1*x1 + B1*C2*x2
                %    [A1 0; 0 B1*C2] * [x1;x2]
                Amerged(stateBinds,feedStateBinds) = ...
                    Amerged(stateBinds,feedStateBinds) ...
                    + B(:,j)*C(outputFeedComp,:);
                
                % part with offset vector k: c1 + B1*k2
                cMerged(stateBinds) = cMerged(stateBinds) ...
                    + B(:,j)*feedFlow.k(outputFeedComp);
                
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

                % sensor noise F: x1'(i) = ... + B1(i,:)*F2(k,:)*v2
                Fmerged(stateBinds,feedNoiseBinds) = ...
                    B(:,j)*feedFlow.F(outputFeedComp,:);
            end
        end
    end
    
    % allocate merged name (remove all default names '')
    namemerged = strjoin(name(cellfun(@(x)~isempty(x),name,'UniformOutput',true)),' x ');
    
    % construct resulting linear systems object
    sys_merged = linearSys(namemerged,Amerged,Bmerged,cMerged,1,[],[],Emerged,Fmerged);
end

function sys_merged = aux_mergeFlowsNonlinearSys(pHA,flowList,locID)
    % note: disturbance/noise matrices not supported!

    % number of components in parallel hybrid automaton
    numComps = length(flowList);
    
    % construct symbolic state vector and input vector for merged flow
    x = sym('x',[pHA.nrOfStates,1]);
    u = sym('u',[pHA.nrOfInputs,1]);
    
    % initialize dynamic function
    f = sym(zeros(pHA.nrOfStates,1));
    
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
    sys_merged = nonlinearSys(name,funHan,pHA.nrOfStates,pHA.nrOfInputs);
end

% ------------------------------ END OF CODE ------------------------------
