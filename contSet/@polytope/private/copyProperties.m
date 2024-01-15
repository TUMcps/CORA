function P_out = copyProperties(P,P_out,mode)
% copyProperties - copies all properties from P to P_out
%
% Syntax:
%    P_out = copyProperties(P,P_out)
%
% Inputs:
%    P - polytope object
%    P_out - polytope object
%    mode - 'all' (copy all properties)
%           'noV' (vertex representation not copied)
%
% Outputs:
%    P_out - polytope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% boundedness
if ~isempty(P.bounded.val)
    P_out.bounded.val = P.bounded.val;
end

% degeneracy
if ~isempty(P.fullDim.val)
    P_out.fullDim.val = P.fullDim.val;
end

% emptiness
if ~isempty(P.emptySet.val)
    P_out.emptySet.val = P.emptySet.val;
end

% minimal representations
if ~isempty(P.minHRep.val)
    P_out.minHRep.val = P.minHRep.val;
end
if ~isempty(P.minVRep.val)
    P_out.minVRep.val = P.minVRep.val;
end

% vertex representation (only if desired)
if ~strcmp(mode,'noV') && ~isempty(P.V.val)
    P_out.V.val = P.V.val;
end

% ------------------------------ END OF CODE ------------------------------
