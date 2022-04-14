function res = mergeTransitionSets(obj, transSets,loc, labelOccs)
% mergeTransitionSets - Compute the transition set of a combined location 
%                       from the transition sets of the subcomponents.
%
%
%
% Syntax:  
%    res = mergeTransitionSets(obj, transSets, loc, inv)
%
% Input:
%     obj - The containing parallelHybridAutomaton
%     transSets - A cell array containing the transition sets for all
%                 components
%     loc - id of the current location 
%     labelOccs - map denoting in how many components a synchronization
%                 label occurs
%
% Outputs:
%    res - cell array containing the resulting transition sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schoepfer, Niklas Kochdumper, Maximilian Perschl
% Written:      14-June-2018
% Last update:  16-March-2022 finalize synchronization implementation
% Last revision: ---

%------------- BEGIN CODE --------------

    numComp = length(transSets);
    res = {};
    cnt = 1;
    % list of checked transition labels
    checkedLabels = [""];

    % loop over all subcomponents
    for i = 1:numComp
        
        compTrans = transSets{i};
        stateBind = obj.bindsStates{i};

        % loop over all transitions for the current subcomponent
        for t = 1:length(compTrans)
            
            trans = compTrans{t};
            
            % prepare target location vector
            targ = loc;
                
            % handle non-synchronized transition
            if isempty(trans.synchLabel)
            
                % update the destination (target idx only changes in component i)
                targ(i) = trans.target;
                
                % convert each component transition to a transition on the PHA
                
                % transitions whose reset depends on input have to be
                % handled seperately for now
                if isfield(trans.reset,'hasInput') && trans.reset.hasInput
                    res{cnt} = projectInputDependentTrans(obj,trans,obj.numStates,i,loc,targ);
                else
                    res{cnt} = projectHighDim(trans,obj.numStates,stateBind,targ);
                end
                cnt = cnt + 1;
            
            % handle synchronized transition
            else
                label = trans.synchLabel;
                % skip over labels that have already been synchronized
                if ismember(label,checkedLabels)
                    continue;
                else
                    checkedLabels = vertcat(checkedLabels,label);
                end
                
                % update target location of component i
                targ(i) = trans.target;

                % get transitions of other components with the same synchronization label
                % start with index i+1, since all labels occuring in components of
                % index <= i have been checked already
                labelTransSet = [trans];
                % along with relevant transitions, we safe from which
                % component they stem to infer the correct state indices
                % later
                labelStateIndices = {};
                labelStateIndices{1} = obj.bindsStates{i};
                for j = i+1:length(transSets)
                    currentComp = transSets{j};
                    for k = 1:length(currentComp)
                        currentTrans = currentComp{k};
                        if strcmp(currentTrans.synchLabel,label)
                            labelTransSet = vertcat(labelTransSet,currentTrans);
                            labelStateIndices = vertcat(labelStateIndices,obj.bindsStates(j));
                            % update target location of component k
                            targ(k) = currentTrans.target;
                            % once a transition with the relevant label has
                            % been found, we can stop scanning the current
                            % component
                            break;
                        end
                    end
                end
                
                % the transition is only active if ALL transitions with the
                % same label across components are active, therefore we
                % check if the number of transitions with the current label
                % is equal to the number of components the label occurs in
                % (if the equality doesn't hold, a transition is blocked
                % because a transition with the same label is in a
                % non-active location of a subcomponent)
                if length(labelTransSet) ~= labelOccs(label)
                    continue;
                end
                
                
                % create a new transition from the transitions to be
                % synchronized, which models the relevant transitions being
                % executed simultaneously 
                
                
                % first, we project each transition to the higher dimension
                for j = 1:length(labelTransSet)
                    if isfield(trans.reset,'hasInput') && trans.reset.hasInput
                        labelTransSet(j) = projectInputDependentTrans(obj,labelTransSet(j),obj.numStates,i,loc,targ);
                    else
                        labelTransSet(j) = projectHighDim(labelTransSet(j),...
                            obj.numStates,labelStateIndices{j},targ);
                    end
                end
                
                % synchronize guards by intersection
                % note: might not work for all nonlinear guards
                
                resultingGuard = labelTransSet(1).guard;
                for j = 2:length(labelTransSet)
                    if ~isempty(labelTransSet(j).guard)
                        if isempty(resultingGuard)
                            resultingGuard = labelTransSet(j).guard;
                        else
                            resultingGuard = resultingGuard & labelTransSet(j).guard;
                        end
                    end
                end
                % reset: merge reset functions (variables are disjoint)
                resultingReset = synchronizeResets(labelTransSet,labelStateIndices,obj.numStates,obj.numInputs);
                
                % add resulting transition to the result
                res{cnt} = transition(resultingGuard,resultingReset,targ);
                cnt = cnt + 1;
            end            
        end
    end
end


% Auxiliary Functions -----------------------------------------------------
function synchedReset = synchronizeResets(transitionSet,stateBindSet,dims,inputDims)
% since the variables of the different subcomponents are considered disjoint in this
% implementation, we can just create a snychronized reset function by
% performing all subcomponent resets in sequence
% (in all reset functions, for the dimensions of states of other components are
%  just the identity function)


    linearResets = zeros(size(transitionSet));
    for i = 1:length(transitionSet)
        if isfield(transitionSet(i).reset,'A')
            linearResets(i) = 1;
        end
    end

    % synchronize purely linear resets
    if all(linearResets)
        A = eye(dims);
        B = zeros(dims,inputDims);
        c = zeros(dims,1);
        % since all resets have been projected to higher dimension, we can
        % just add them
        for i = 1:length(transitionSet)
            A = A + transitionSet(i).reset.A;
            if isfield(transitionSet(i).reset,'B')
                B = B + transitionSet(i).reset.B;
            end
            c = c + transitionSet(i).reset.c;
        end
        synchedReset.A = A;
        synchedReset.B = B;
        synchedReset.c = c;
        return;
        
    % convert linear resets to nonlinear resets
    else
        for i = 1:length(linearResets)
            if linearResets(i)
                if isfield(transitionSet(i).reset,'B')
                    new_f = @(x,u) transitionSet(i).reset.A * x + transitionSet(i).reset.B * u + transitionSet(i).reset.c;
                    newReset.hasInput = 1;
                    newReset.inputDim = size(transitionSet(i).reset.B,2);

                else
                    new_f = @(x) transitionSet(i).reset.A * x + transitionSet(i).reset.c;
                    newReset.hasInput = 0;
                end
                newReset.f = new_f;
                newTransition = transition(transitionSet(i).guard,newReset,transitionSet(i).target, ...
                    transitionSet(i).synchLabel,size(transitionSet(i).reset.A,1));
                transitionSet(i) = newTransition;
            end
        end
    end

    % synchronize nonlinear resets
    % procedure: merge the resets and jacobians, 
    % then build the sums of hessians and ternary tensors
        
    % nonlinear:
    synchedReset.f = @(x) x;
    for j = 1:length(transitionSet)
        synchedReset.f = @(x) transitionSet(j).reset.f(synchedReset.f(x));
    end
    
    synchedReset.J = @(x) mergeJacobians(x,obj.numStates,transitionSet,stateBindSet);
    
    temp = @(x) zeros(dims);
    synchedReset.Q = repmat({temp},[dims,1]);
    
    % Hessians
    for i = 1:dims
        for j = 1:length(transitionSet)
            synchedReset.Q{i} = @(x)synchedReset.Q{i}(x) + transitionSet(j).reset.Q{i}(x);
        end
    end
    
    
    synchedReset.T = repmat({temp},[dims,dims]);
    % Third order tensors
    for i = 1:dims
        for j = 1:dims
            for k = 1:length(transitionSet)
                if ~isempty(transitionSet(k).reset.T{i,j})
                    synchedReset.T{i}{j} = @(x)synchedReset.T{i,j}(x) + transitionSet(k).reset.T{i,j}(x);
                end
            end
        end
    end   
end


function res = mergeJacobians(x,n,transitionSet,stateBindSet)
    res = eye(n); 
    for i = 1:length(transitionSet)
        res(stateBindSet(i),stateBindSet(i)) = transitionSet(i).reset.J(x(stateBindSet(i)));
    end
end

% function res = mergeQuadMat(x,n,transitionSet,stateBindSet)
%     res = zeros(n);
%     for i = 1:length(transitionSet)
%         res(stateBindSet(i),stateBindSet(i)) = transitionSet(i).Q(x(stateBindSet(i)));
%     end
% end
% 
% function res = mergeThirdOder(x,n,transitionSet,stateBindSet)
%     res = interval(zeros(n));
%     for i = 1:length(transitionSet)
%         res(stateBindSet(i),stateBindSet(i)) = transitionSet(i).T(x(stateBindSet(i)));
%     end
% end
    

%------------- END OF CODE --------------