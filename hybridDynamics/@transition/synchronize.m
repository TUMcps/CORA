function trans = synchronize(transList,idStates,locID,compIdx,stateBinds,inputBinds,flowList)
% synchronize - synchronizes multiple transitions: the guard sets are
%    intersected, the reset functions are combined into a single function,
%    and the individual targets are combined to a joint updated target
%
% Syntax:
%    trans = synchronize(transList,idStates,locID,compIdx,stateBinds,inputBinds,flowList)
%
% Inputs:
%    transList - class array of transition objects
%    idStates - states whose reset is mapped by identity
%    locID - IDs of the currently active locations
%    compIdx - components corresponding to elements in transList
%    stateBinds - states of the high-dimensional space that correspond to
%                 the states of the low-dimensional reset object
%    inputBinds - connections of input to global input/outputs of other
%                 components
%    flowList - list of flow equations for each component
%
% Outputs:
%    trans - resulting transition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that...
% 1. guard sets are of same dimension
% 2. reset functions are of same dimension
% 3. target vectors are of equal size

% synchronize guards by intersection
guard = fullspace(dim(transList(1).guard));
for i=1:length(transList)
    guard = and_(guard,transList(i).guard,'exact');
end

% synchronize resets (without resolution of inputs), but they must all be
% of the same class, so we may need to convert before synchronization
resets = aux_convertResets(transList);
reset_sync = synchronize(resets,idStates);

% note: number of total cells in outputEq must equal reset.inputDim
reset = resolve(reset_sync,flowList,stateBinds,inputBinds);

% compose the target vector by overriding all entries in the currently
% location ID vector for which there is a transition
target = locID;
for i=1:length(compIdx)
    target(compIdx) = transList(i).target;
end

% instantiate resulting transition (no synchronization label anymore)
trans = transition(guard,reset,target);

end


% Auxiliary functions -----------------------------------------------------

function resets = aux_convertResets(transList)

% check which ones are linearReset and which are nonlinearReset objects 
isLinearReset = arrayfun(@(trans) isa(trans.reset,'linearReset'),...
    transList,'UniformOutput',true);
isNonlinearReset = arrayfun(@(trans) isa(trans.reset,'nonlinearReset'),...
    transList,'UniformOutput',true);

if all(isLinearReset) || all(isNonlinearReset)
    resets = [transList.reset];
    return
end

% convert all linear reset functions to nonlinear reset functions: since
% the order of reset functions does not matter (they must have been
% correctly lifted to the same dimension with correct binding before), we
% first gather all nonlinear reset functions and then append the remaining
% ones (note: isNonlinearReset has at least one true element)
resets = [transList(isNonlinearReset).reset];
resets_converted = arrayfun(@(r) nonlinearReset(r),...
    transList(isLinearReset).reset,'UniformOutput',true);
resets = [resets; resets_converted];

end

% ------------------------------ END OF CODE ------------------------------
