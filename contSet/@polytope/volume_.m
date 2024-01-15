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

% ------------------------------ BEGIN CODE -------------------------------

% fully empty
if representsa_(P,'fullspace',0)
    vol = Inf; return
elseif representsa(P,'emptySet',1e-10)
    vol = 0; return
end

% 1D case very fast
if dim(P) == 1
    if ~isempty(P.V.val)
        V = P.V.val; 
    else
        V = vertices(P);
    end
    
    vol = max(V) - min(V);
    return
end

%compute vertices
try
    if ~isempty(P.V.val)
        V = P.V.val; 
    else
        V = vertices(P);
        P.V.val = V;
    end
catch ME
    vol = aux_specialCases(P,ME);
    return  
end

if size(V,2) <= dim(P)
    % a polytope needs at least n+1 vertices to be full-dimensional
    vol = 0;

elseif size(V, 1) == dim(P)+1
    % cheap calculation for full-dimensional simplex
    S = V';
    D = zeros(size(S, 1), size(S, 2)-1);
    for j = 2:size(S, 2)
		D(:, j-1) = S(:, j)-S(:, 1);
    end
    vol = 1/factorial(size(S, 2)-1)*abs(det(D));

else
    % general computation
    try
        [~, vol] = convhulln(V');
    catch ME
        vol = aux_specialCases(P,ME);
        return
    end

end

end


% Auxiliary functions -----------------------------------------------------

function vol = aux_specialCases(P,ME)

    % check if empty or not full dimensional
    if representsa(P, 'emptySet') || ~isFullDim(P)
        vol = 0;
        return
    end
    
    %check if unbounded
    if ~isBounded(P)
        vol = Inf;
        return
    end

    % nothing we can do
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
