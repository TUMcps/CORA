function list = indexList(loc)
% indexList - returns a list that relates active events to guards; an
%    invariant is defined as 0, other numbers refer to the guard number
%
% Syntax:
%    list = indexList(loc)
%
% Inputs:
%    loc - location object
%
% Outputs:
%    list - list of guard numbers, position refers to the event
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-May-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % get indices of the invariant
    eq = aux_getEquations(loc.invariant);
    list(1:eq) = 0;

    % get indices of the guards
    for i = 1:length(loc.transition)
        eq = aux_getEquations(loc.transition(i).guard);
        list((end+1):(end+eq)) = i;
    end
end


% Auxiliary functions -----------------------------------------------------

function eq = aux_getEquations(S)
% get number of inequality constraints describing the set representation

    if isa(S,'polytope')

        eq = length(S.b);
    
    elseif isa(S,'levelSet')
    
        if iscell(S.compOp)
            eq = length(S.compOp);
        else
            eq = 1; 
        end
    
    else
        throw(CORAerror('CORA:specialError',...
            'Guard set has to be polytope object or levelSet object'));
    end
end

% ------------------------------ END OF CODE ------------------------------
