function invariant_cora2spaceex(Obj,location, docNode, inv)
% invariant_cora2spaceex - 
%
% Syntax:
%    invariant_cora2spaceex(Obj,location, docNode, inv)
%
% Inputs:
%    Obj -
%    location -
%    docNode -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cora2spaceex

% Authors:       Farah Atour
% Written:       24-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isa(Obj,'hybridAutomaton')
    
        % convert invariant set to polytope if possible
    if isa(inv,'interval') || isa(inv,'zonotope') || ...
       isa(inv,'conZonotope') || isa(inv,'zonoBundle')
       inv = polytope(inv);
    end
    
    % convert the different types of invariant sets
    if isa(inv,'conHyperplane')
        
        eqs = conHyperplane_cora2spaceex(inv);
        
    elseif isa(inv,'polytope')
        
        eqs = polytope_cora2spaceex(inv);
        
    elseif isa(inv,'levelSet')
        
        eqs = levelSet_cora2spaceex(inv);
        
    else     
        throw(CORAerror('CORA:notSupported','Given invariant set class not supported.'));
    end
    
    % Add the element node (invariant), for the parent element (location) and
    % set the equation attribute.
    invariant = docNode.createElement('invariant');
    invariant.appendChild(docNode.createTextNode(eqs));
    location.appendChild(invariant);
    
else
    
    %Add the empty element node (invariant), for the parent element (location)
    invariant = docNode.createElement('invariant');
    location.appendChild(invariant);
    
end

end

% ------------------------------ END OF CODE ------------------------------
