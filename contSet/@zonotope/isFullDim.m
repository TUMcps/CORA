function [res, subspace] = isFullDim(Z,varargin)
% isFullDim - checks if the dimension of the affine hull of a zonotope is
%    equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(Z)
%    res = isFullDim(Z,tol)
%    [res,subspace] = isFullDim(_)
%
% Inputs:
%    Z - zonotope object
%    tol - numeric, tolerance
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that Z is strictly contained in
%               center(Z)+span(x_1,...,x_k)
%               (here, 'strictly' means that k is minimal).
%               Note that if Z is just a point, subspace=[].
%
% Example:
%    Z1 = zonotope([1 2 1;3 1 2]);
%    Z2 = zonotope([1 2 1;3 4 2]);
%
%    isFullDim(Z1)
%    isFullDim(Z2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Adrian Kulmburg
% Written:       02-January-2020 
% Last update:   12-March-2021 (MW, add empty case)
%                04-February-2024 (AK, add subspace computation)
%                17-May-2024 (TL, added tol)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
tol = setDefaultValues({1e-6},varargin);

if ~representsa_(Z,'emptySet',eps)
    % If the user only wants to know if the set is non-degenerate, we can do
    % this quickly with the rank (note that this requires the computation
    % of the SVD)
    if nargout < 2
        Zdim = dim(Z);
        Grank = rank(Z.G,tol);
        res = Zdim == Grank;
    else
        % If the user also wants the subspace, we need to make a SVD no matter
        % what, so might as well check for the rank ourselves
        Zdim = dim(Z);
    
        [U,Sigma, ~] = svd(Z.G);
    
        s = diag(Sigma);
        r = sum(s>tol);
    
        res = Zdim == r;
        if ~res
            subspace = U(:,1:r);
        else
            subspace = eye(Zdim);
        end

    end
else
    res = false;
    subspace = [];
end

% ------------------------------ END OF CODE ------------------------------
