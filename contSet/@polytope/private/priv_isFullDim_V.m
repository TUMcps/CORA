function [res,subspace] = priv_isFullDim_V(V,tol)
% priv_isFullDim_V - checks if the dimension of the affine hull of a 
%    polytope is equal to the dimension of its ambient space; additionally,
%    one can obtain a basis of the subspace in which the polytope is
%    contained
%
% Syntax:
%    res = priv_isFullDim_V(V,tol)
%    [res,subspace] = priv_isFullDim_V(V,tol)
%
% Inputs:
%    V - vertex representation
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that P is strictly contained in
%               center(P)+span(x_1,...,x_k)
%               (here, 'strictly' means that k is minimal).
%               Note that if P is just a point, subspace=[].
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

% read out dimension
n = size(V,1);

% compute rank shifted by mean
V_shifted = V - mean(V,2);
rankV = rank(V_shifted,tol);
% compare to ambient dimension
res = rankV == n;

% compute basis of affine hull
if nargout == 2
    if res
        subspace = eye(n);
    else
        [Q,R] = qr(V_shifted);
        subspace = Q(:,1:rankV);
    end
else
    % dummy value for subspace (not used further)
    subspace = NaN;
end

% ------------------------------ END OF CODE ------------------------------
