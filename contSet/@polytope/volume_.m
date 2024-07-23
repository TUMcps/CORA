function vol = volume_(P,varargin)
% volume_ - computes the volume of a polytope
%
% Syntax:
%    vol = volume_(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    vol - volume
%
% Example: 
%    P = polytope([1 1; -1 0; 0 -1],[1; 0; 0]);
%    vol = volume(P)
%
% References: MPT-toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Viktor Kotsev
% Written:       20-June-2022
% Last update:   ---
% Last revision: 25-May-2023 (MW, speed up by throwing errors later)
%                12-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% cheap shots at the start: unbounded, empty, degenerate
if (~isempty(P.bounded.val) && ~P.bounded.val) || representsa_(P,'fullspace',0)
    vol = Inf; return
elseif ~isempty(P.emptySet.val) && P.emptySet.val
    vol = 0; return
elseif ~isempty(P.fullDim.val) && ~P.fullDim.val
    vol = 0; return
end

% read out dimension
n = dim(P);

% compute vertices (note: they are read out if already computed)
try
    V = vertices_(P,'lcon2vert');
catch ME
    vol = aux_volume_specialCases(P,ME);
    return
end

% different cases
if n == 1
    % 1D: volume is the length of the line between min and max vertices
    vol = max(V) - min(V);

elseif size(V,2) <= n || (~isempty(P.fullDim.val) && ~P.fullDim.val)
    % a polytope needs at least n+1 vertices to be full-dimensional;
    % degeneracy may also have been detected in vertex enumeration above
    vol = 0;

elseif size(V, 1) == n+1
    % cheap calculation for full-dimensional simplex
    S = V';
    D = zeros(size(S,1), size(S,2)-1);
    for j = 2:size(S,2)
		D(:,j-1) = S(:,j)-S(:,1);
    end
    vol = 1/factorial(size(S,2)-1)*abs(det(D));

else
    % general computation via additional output argument of Quickhull
    try
        [~, vol] = convhulln(V');
    catch ME
        vol = aux_volume_specialCases(P,ME);
        return
    end

end

end


% Auxiliary functions -----------------------------------------------------

function vol = aux_volume_specialCases(P,ME)

    % check if empty or not full dimensional
    if representsa(P,'emptySet') || ~isFullDim(P)
        vol = 0;
        return
    end
    
    % check if unbounded
    if ~isBounded(P)
        vol = Inf;
        return
    end

    % nothing we can do
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
