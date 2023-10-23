function guard_cora2spaceex(tran, docNode, guard)
% guard_cora2spaceex - 
%
% Syntax:
%    guard_cora2spaceex(tran, docNode, guard)
%
% Inputs:
%    tran - transition
%    docNode - 
%    guard - guard set
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


    % convert guard sets to polytopes if possible
    if isa(guard,'interval') || isa(guard,'zonotope') || ...
       isa(guard,'conZonotope') || isa(guard,'zonoBundle')
       guard = polytope(guard);
    end
    
    % convert the different types of guard sets
    if isa(guard,'conHyperplane')
        
        eqs = conHyperplane_cora2spaceex(guard);
        
    elseif isa(guard,'polytope')
        
        eqs = polytope_cora2spaceex(guard);
        
    elseif isa(guard,'levelSet')
        
        eqs = levelSet_cora2spaceex(guard);
        
    elseif isnumeric(guard) && isempty(guard)
        % TODO: replace by fullspace/emptySet
        eqs = '';

    else
        throw(CORAerror('CORA:converterIssue',...
            'Given guard set class not supported.'));
    end

    % Add the element node (guard), for the parent element (transition) and
    % set the equation attribute.
    guard = docNode.createElement('guard');
    guard.appendChild(docNode.createTextNode(eqs));
    tran.appendChild(guard);
end

% ------------------------------ END OF CODE ------------------------------
