function V = removeDuplicates(V,varargin)
% removeDuplicates - remove duplicates in an array of vertices up to a
%    given tolerance
%
% Syntax:
%    V = removeDuplicates(V)
%    V = removeDuplicates(V,tol)
%
% Inputs:
%    V - vertices
%    tol - (optional) tolerance
%
% Outputs:
%    V - vertices
%
% Example:
%    V = 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set default tolerance
tol = setDefaultValues({1e-10}, varargin);

% number of vertices
numVert = size(V,2);

% too few vertices
if size(V,2) < 2
    return
end

% indices for which vertices are kept
idxKept = true(1,numVert);

% check distance between vertices
for i = 1:numVert
    % skip if i-th vertex already removed
    if ~idxKept(i)
        continue
    end

    % compute 2-norm of i-th vertex to all others
    dist = vecnorm(V - V(:,i));

    % indices where distance is smaller than tolerance
    idxRemove = dist < tol;
    % eliminate vertices that are already selected for removal
    idxRemove = idxRemove & idxKept;

    if nnz(idxRemove) > 1
        % note: idxRemove is always at least 1 since the distance of the
        % vertex to itself is 0

        % replace i-th vertex by mean of all vertices that are in its
        % neighborhood
        V(:,i) = mean(V(:,idxRemove),2);

        % update indices
        idxRemove(i) = false;
        idxKept = idxKept & ~idxRemove;
    end
end

% remove duplicate vertices (if any)
if ~all(idxKept)
    V = V(:,idxKept);
end

% ------------------------------ END OF CODE ------------------------------
