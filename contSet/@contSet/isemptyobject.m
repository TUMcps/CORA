function res = isemptyobject(S)
% isemptyobject - checks whether a contSet object contains any information
%    at all; consequently, the set is equivalent to the empty set 
%
% Syntax:
%    res = isemptyobject(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    isemptyobject(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-June-2022
% Last update:   18-August-2022 (MW, extend to class arrays)
%                24-July-2023 (MW, move checks to subclasses, throw error)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:noops',S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(S)
% properties must be same as in fully-empty object of same class
% this is a rather brute-force version that needs to be updated
% semi-regularly, but for efficiency reasons there is no real other option

classname = class(S);

% different checks depending on class
switch classname

    case 'Polyhedron'
        % from mpt toolbox
        res = islogical(S.irredundantVRep) && ~S.irredundantVRep ...
            && islogical(S.irredundantHRep) && ~S.irredundantHRep ...
            && islogical(S.hasHRep) && ~S.hasHRep ...
            && islogical(S.hasVRep) && ~S.hasVRep ...
            && isnumeric(S.A) && isempty(S.A) ...
            && isnumeric(S.b) && isempty(S.b) ...
            && isnumeric(S.Ae) && isempty(S.Ae) ...
            && isnumeric(S.be) && isempty(S.be) ...
            && isnumeric(S.H) && isempty(S.H) ...
            && isnumeric(S.He) && isempty(S.He) ...
            && isnumeric(S.R) && isempty(S.R) ...
            && isnumeric(S.V) && isempty(S.V) ...
            && isnumeric(S.Dim) && isscalar(S.Dim) && S.Dim == 0 ...
            && isnumeric(S.Data) && isempty(S.Data);

<<<<<<< HEAD
    case 'polytope'
        res = isnumeric(S.A) && isempty(S.A) ...
            && isnumeric(S.b) && isempty(S.b) ...
            && isnumeric(S.Ae) && isempty(S.Ae) ...
            && isnumeric(S.be) && isempty(S.be);

 
end

end

function res = aux_checkIfEmpty_old(S,S_empty,props)
% old method based on instantiated empty object (which takes a long time,
% therefore replaced by more direct albeit less sustainable solution)
% this function is kept for debugging purposes;
% note: compute other two input arguments by
%    S_empty = eval([classname '();']);
%    props = properties(S);

% properties must be same as in fully-empty object of same class

    % init return value
    res = true;

    nrProps = length(props);
    for p=1:nrProps
    
        % property of contSet object
        prop = S.(props{p});
    
        % corresponding property of fully-empty set
        prop_empty = S_empty.(props{p});
    
        % type of property
        type = class(prop);
    
        % check depends on type
        switch type
    
            case 'logical'
                if prop ~= prop_empty
                    res = false; return
                end
    
            case {'double','cell','sym'}
                if ~isequal(prop,prop_empty)
                    res = false; return
                end
    
            case {'char','string'}
                if ~strcmp(prop,prop_empty)
                    res = false; return
                end
    
            otherwise
                % check if the property is a contSet object (e.g., halfspace
                % object in conHyperplane object)
                if isa(prop,'contSet')
                    % recursive call
                    if ~isemptyobject(prop)
                        res = false; return
                    end
                else
                    throw(CORAerror('CORA:specialError','Unknown data type'));
                end
        end
        
    end

end

% ------------------------------ END OF CODE ------------------------------
