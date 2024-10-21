function V = priv_compact_V(V,tol)
% priv_compact_V - removes redundancies in the vertex representation
%
% Syntax:
%    V = priv_compact_V(V,tol)
%
% Inputs:
%    V - vertex representation
%    tol - tolerance
%
% Outputs:
%    V - vertex representation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract arguments
[n,numVert] = size(V);

% if at most one vertex, nothing to do...
if numVert <= 1
    return;
end

% handle 1-dimensional case separately
if n == 1
    min_V = min(min(V));
    max_V = max(max(V));
    V = [min_V; max_V]';    
    return
end

% compute convex hull and take only the unique points
try
    % If the polytope is degenerate, we first need to restrict it
    % to its affine hull
    [ifd, B] = priv_isFullDim_V(V,tol);
    if ifd
        K = convhulln(V');
        V = V(:,unique(K));
    else
        % center of the polytope
        c = sum(V,2)./size(V,2);

        % decompose and compute the rank of B
        [U, ~, ~] = svd(B);
        r = size(B,2);

        V_affSubspace = [eye(r) zeros([r (n-r)])] * U' * (V - c);
        if size(V_affSubspace,1) == 1
            % convhulln does not work for 1D sets of vertices...
            minV = min(V_affSubspace);
            maxV = max(V_affSubspace);
            if ~withinTol(minV,maxV,tol)
                V_affSubspace = [minV, maxV];
            else
                V_affSubspace = minV;
            end
        else
            K_affSubspace = convhulln(V_affSubspace');
            V_affSubspace = V_affSubspace(:, unique(K_affSubspace));
        end

        % back-projection
        V = U * [eye(r); zeros([(n-r) r])] * V_affSubspace + c;
    end
catch ME
    if strcmp(ME.identifier,'MATLAB:cgprechecks:NotEnoughPts')
        % not enough unique points specified -> we assume that this
        % means minimal representation
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
