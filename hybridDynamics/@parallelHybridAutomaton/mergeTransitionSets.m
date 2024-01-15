function [transSet,mergedLabels] = mergeTransitionSets(pHA,transList,locID,allLabels)
% mergeTransitionSets - computes the transition set of a combined location 
%    from the transition sets of the subcomponents
%
% Syntax:
%    transSet = mergeTransitionSets(pHA,transList,locID,allLabels)
%    [transSet,mergedLabels] = mergeTransitionSets(pHA,transList,locID,allLabels)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%    transList - A cell array containing the transition sets for all
%                 components
%    locID - IDs of the currently active locations
%    allLabels - information about all synchronization labels
%
% Outputs:
%    transSet - cell array containing the resulting transition sets
%    mergedLabels - list of synchronization labels which were used to
%                   construct a merged transition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton/mergeInvariants

% Authors:       Johann Schoepfer, Niklas Kochdumper, Maximilian Perschl
% Written:       14-June-2018
% Last update:   16-March-2022 (MP, implement synchronization labels)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % read out number of components
    numComp = length(pHA.components);

    % initialize resulting set of transitions
    transSet = transition();
    % counter for number of transitions (cannot be known beforehand due to
    % synchronization labels)
    cnt = 1;
    
    % list of checked transition labels
    checkedLabels = {};
    % list of synchronization labels that were merged
    mergedLabels = {};

    % loop over all subcomponents
    for i = 1:numComp
        
        compTrans = transList{i};
        stateBind = pHA.bindsStates{i};

        % loop over all transitions for the current subcomponent
        for t = 1:length(compTrans)
            
            trans = compTrans(t);
            
            % prepare target location vector
            resultingTarget = locID;
                
            % handle non-synchronized transition (no synchronization label)
            if isempty(trans.syncLabel)
            
                % update the destination (target idx only changes in component i)
                resultingTarget(i) = trans.target;
                
                % convert each component transition to a transition on the PHA
                
                % transitions whose reset depends on input have to be
                % handled seperately for now; here, the other states in the
                % projected reset function are set to identity
                if trans.reset.hasInput
                    transSet(cnt) = projectInputDependentTrans(pHA,...
                        trans,pHA.dim,i,locID,resultingTarget,true);
                else
                    transSet(cnt) = lift(trans,pHA.dim,...
                        stateBind,resultingTarget,true);
                end

                % increment counter of resulting transitions
                cnt = cnt + 1;
            
            % handle synchronized transition
            else
                syncLabel = trans.syncLabel;
                % skip over labels that have already been synchronized;
                % additionally, we can skip labels where there is not one
                % active instance per component containing that label
                if ismember(syncLabel,checkedLabels)
                    continue;
                else
                    checkedLabels = [checkedLabels; syncLabel];
                end

                % for partly-synchronized transitions, the states of all
                % components, which do not participate in the
                % synchronization, remain (identity)
                idStates = zeros(pHA.dim,1);
                
                % update target location of component i
                resultingTarget(i) = trans.target;

                % get transitions of other components with the same label:
                % start with index i+1 as all labels occurring in
                % components of index <= i have been checked already
                labelTransSet = trans;
                labelCompIdx = i;

                % we also save from which component they originate to infer
                % the correct state indices later
                labelStateIndices = {};
                labelStateIndices{1} = pHA.bindsStates{i};

                % loop over each component: cannot start at component i+1
                % (which would suffice because otherwise the label in the
                % current component would have been checked before), since
                % we need the states with identity in the reset
                for j = 1:length(transList)
                    if j ~= i
                        idx = arrayfun(@(x) strcmp(x.syncLabel,syncLabel),...
                            transList{j},'UniformOutput',true);
                        if ~isempty(idx) && any(idx)
                            % only one transition per component may have a
                            % given synchronization label, therefore we can use
                            % this index to read out the transition and other
                            % relevant information
                            temp = transList{j}(idx);
                            labelTransSet = [labelTransSet; temp(1)];
                            labelCompIdx = [labelCompIdx; j];
                            labelStateIndices = [labelStateIndices; pHA.bindsStates(j)];
    
                            % update target location of component j
                            resultingTarget(j) = temp(1).target;
                        else
                            % set indices of state binds to 1 (will later be
                            % used in reset function to maintain values)
                            idStates(pHA.bindsStates{j}) = 1;
                        end
                    end
                end
                
                % the transition is only active if in all components where
                % the label occurs, one of these transitions is active
                % (the lenghty variable on the right-hand side computes the
                % number of components in which <label> occurs)
                if length(labelTransSet) ~= ...
                        length(unique(allLabels(strcmp({allLabels.name},syncLabel)).component))
                    continue;
                end
                
                
                % create a new transition from the transitions to be
                % synchronized, which models the relevant transitions being
                % executed simultaneously

                % save label for product of invariants
                mergedLabels = [mergedLabels; syncLabel];

                
                % first, we project each transition to the full dimension;
                % here, the other states of the reset function are set to
                % zero since they will be either handled by their own reset
                % function or if there is none, then in synchronizeResets
                for j = 1:length(labelTransSet)
                    if labelTransSet(j).reset.hasInput
                        labelTransSet(j) = projectInputDependentTrans(pHA,...
                            labelTransSet(j),pHA.dim,labelCompIdx(j),...
                            locID,resultingTarget,false);
                    else
                        labelTransSet(j) = lift(labelTransSet(j),...
                            pHA.dim,labelStateIndices{j},resultingTarget,false);
                    end
                end
                % for nonlinear resets: all resulting reset functions are
                % defined using the full state as reset.stateDim, the
                % number of global inputs as reset.inputDim (so that
                % reset.hasInput is only true if there are global inputs to
                % composed system), and reset.f takes the extended state
                % vector (full state + global inputs) as input argument
                % for linear resets: states use the A matrix (where inputs
                % which are states from other components are already
                % resolved), global inputs use the B matrix
                
                % synchronize guards by intersection
                % note: might not work for all nonlinear guards
                resultingGuard = fullspace(pHA.dim);

                % note: empty guard sets can simply be skipped as they
                % represent the full subspace of the covered dimensions,
                % which is conveniently represented by the resulting zeros
                % -> the guard intersection will correctly yield
                % full-dimensional sets along these dimensions, thus we do
                % not require any further processing at this stage
                for j=1:length(labelTransSet)
                    if ~isa(labelTransSet(j).guard,'fullspace')
                        if isa(resultingGuard,'fullspace')
                            % first full-dimensional guard enters here
                            resultingGuard = labelTransSet(j).guard;
                        else
                            % all subsequent full-dimensional guards
                            resultingGuard = and_(resultingGuard,labelTransSet(j).guard,'exact');
                        end
                    end
                end
                

                % merge reset functions into one reset function
                resultingReset = synchronizeResets(labelTransSet,...
                    pHA.dim,pHA.nrOfInputs,idStates);
                
                % instantiate resulting transition
                transSet(cnt) = transition(resultingGuard,resultingReset,...
                    resultingTarget);

                % increment counter for number of transitions
                cnt = cnt + 1;

            end            
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
