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
%    ap - atomicProposition object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

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
    function ap = atomicProposition(set, dims, locs)
        % constructor 
        ap.set = set;

        if nargin < 3
            ap.locs = [];
        else
            ap.locs = locs;
        end

        if nargin < 2
            ap.dims = 1:dim(set);
        else
            ap.dims = dims;
        end
    end

    function out = containsLoc(ap, loc)
        % check if loc is contained
        if ~isempty(ap.locs)
            out = in(ap.locs, loc);
        else
            out = true;
        end
    end

    function out = evaluatePoint(ap, point, loc)
        % Check if the set contains the projection of the given point.
        out = containsLoc(ap, loc) && contains(ap.set, point(ap.dims));
    end

    function out = canBeTrue(ap, set, loc)
        % check if loc is contained at all
        if ~containsLoc(ap, loc)
            out = false; return
        end

        % If the set is intersecting with the proposition, it can be true.
        set = project(set, ap.dims);
        tol = 1e-12;
        try
            out = isIntersecting_(ap.set, set, 'exact', tol);
        catch e
            if ismember(e.identifier, {'CORA:noops','CORA:noExactAlg'})
                out = isIntersecting_(ap.set, set, 'approx', tol);
            else
                rethrow(e);
            end
        end
    end

    function out = canBeFalse(ap, set, loc)
        % If the set is not enclosed by the proposition, it can be false.
        set = project(set, ap.dims);

        out = ~containsLoc(ap, loc);

        try
            out = out || ~contains(ap.set, set);
        catch e
            if ismember(e.identifier, {'CORA:noops','CORA:noExactAlg'})
                out = out || ~contains(ap.set, set, 'approx');
            else
                rethrow(e);
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
