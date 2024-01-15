classdef atomicProposition
% atomicProposition - geometric representation of an atomic proposition
%
% Syntax:
%    ap = atomicProposition(set)
%    ap = atomicProposition(set, dims)
%    ap = atomicProposition(set, dims, locs)
%
% Inputs:
%    set - geometric representation of the set
%    dims - list of dimensions for this proposition
%    locs - list of locations that are allowed (hybrid automata)
%
% Outputs:
%    obj - atomicProposition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Benedikt Seidl
% Written:       07-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % geometric representation
    set

    % list of dimensions
    dims {mustBePositive}

    % list of locations
    locs {mustBePositive}
end

methods
    function obj = atomicProposition(set, dims, locs)
        obj.set = set;

        if nargin < 3
            obj.locs = [];
        else
            obj.locs = locs;
        end

        if nargin < 2
            obj.dims = 1:dim(set);
        else
            obj.dims = dims;
        end
    end

    function out = containsLoc(obj, loc)
        if ~isempty(obj.locs)
            out = in(obj.locs, loc);
        else
            out = true;
        end
    end

    function out = evaluatePoint(obj, point, loc)
        % Check if the set contains the projection of the given point.

        out = containsLoc(obj, loc) && contains(obj.set, point(obj.dims));
    end

    function out = canBeTrue(obj, set, loc)
        % If the set is intersecting with the proposition, it can be true.
        set = project(set, obj.dims);

        out = containsLoc(obj, loc);

        try
            out = out && isIntersecting(obj.set, set);
        catch e
            if ismember(e.identifier, {'CORA:noops','CORA:noExactAlg'})
                out = out && isIntersecting(obj.set, set, 'approx');
            else
                rethrow(e);
            end
        end
    end

    function out = canBeFalse(obj, set, loc)
        % If the set is not enclosed by the proposition, it can be false.
        set = project(set, obj.dims);

        out = ~containsLoc(obj, loc);

        try
            out = out || ~contains(obj.set, set);
        catch e
            if ismember(e.identifier, {'CORA:noops','CORA:noExactAlg'})
                out = out || ~contains(obj.set, set, 'approx');
            else
                rethrow(e);
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
