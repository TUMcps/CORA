function res = mergeTransitionSets(obj, transSets,loc)
% mergeTransitionSets - Compute the transition set of a combined location 
%                       from the transition sets of the subcomponents.
%
% This composition function does NOT perform any kind of synchronization.
% (transiton.input/ouputLabel fields are unused as of this writing)
%
% Syntax:  
%    res = mergeTransitionSets(obj, transSets, loc, inv)
%
% Input:
%     obj - The containing parallelHybridAutomaton
%     transSets - A cell array containing the transition sets for all
%                 transitions
%     loc - id of the current location 
%
% Outputs:
%    res - cell array containing the resulting transition sets
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Johann Schöpfer, Niklas Kochdumper
% Written:      14-June-2018
% Last update:  09-July-2018 (NK, use function "projectHighDim")
% Last revision: ---

%------------- BEGIN CODE --------------

    numComp = length(transSets);
    res = {};
    counter = 0;
    
    % loop over all subcomponents
    for i = 1:numComp
        
        compTrans = transSets{i};
        stateBind = obj.bindsStates{i};

        % loop over all transitions for the current subcomponent
        for t = 1:length(compTrans)
            
            % convert each component transition to a transition on the PHA
            trans = compTrans{t};

            % project guard set to the higher dimension
            guard = trans.guard;
            if ~isa(guard,'levelSet') && ~isa(guard,'mptPolytope') && ...
               ~isa(guard,'conHyperplane')
                guard = mptPolytope(guard);
            end
            
            guardProj = projectHighDim(guard,obj.numStates,stateBind);

            % project reset function to the higher dimension
            resetProj = projectReset(trans.reset,stateBind,obj.numStates);

            % update the destination (target idx only changes in component i)
            targetProj = loc;
            targetProj(i) = trans.target;

            % construct transition object & add it to set
            counter = counter + 1;
            res{counter} =...
                transition(guardProj,resetProj,targetProj);
        end
    end
end

%------------- END OF CODE --------------