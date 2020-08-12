function [list] = indexList(obj)
% indexList - returns a list that relates active events to guards; an
% invariant is defined as 0, other numbers refer to the guard number
%
% Syntax:  
%    [list] = indexList(obj)
%
% Inputs:
%    obj - location object
%
% Outputs:
%   list - list of guard numbers, whereas the list position refers to the
%   event
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 07-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % get indices of the invariant
    eq = getEquations(obj.invariant);
    list(1:eq) = 0;

    % get indices of the guards
    for i = 1:length(obj.transition)
        eq = getEquations(obj.transition{i}.guard);
        list((end+1):(end+eq)) = i;
    end
end


% Auxiliary Functions -----------------------------------------------------

function eq = getEquations(set)
% get number of inequality constraints describing the set representation

    if isa(set,'mptPolytope')
       
        eq = length(get(set,'b'));
        
    elseif isa(set,'levelSet')
       
       if iscell(set.compOp)
          eq = length(set.compOp);
       else
          eq = 1; 
       end
        
    else
        error('Something went wrong!');      
    end
end

%------------- END OF CODE --------------