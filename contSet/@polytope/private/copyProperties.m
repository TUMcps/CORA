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

% H representation given
if ~isempty(P.isHRep.val)
    P_out.isHRep.val = P.isHRep.val;
end

% V representation given
if ~isempty(P.isVRep.val)
    P_out.isVRep.val = P.isVRep.val;
end

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
if ~strcmp(mode,'noV') && P.isVRep
    P_out.V = P.V;
end

% ------------------------------ END OF CODE ------------------------------
