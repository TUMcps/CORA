function [res,subspace] = isFullDim(cZ,varargin)
% isFullDim - checks if the dimension of the affine hull of a constrained
%    zonotope is equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(cZ)
%    res = isFullDim(cZ,tol)
%    [res,subspace] = isFullDim(cZ)
%
% Inputs:
%    cZ - conZonotope object
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    subspace - (optional) Returns a set of orthogonal unit vectors
%               x_1,...,x_k such that cZ is strictly contained in
%               center(cZ)+span(x_1,...,x_k)
%               (here, 'strictly' means that k is minimal).
%               Note that if cZ is just a point, subspace=[].
%
% Example:
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%
%    P = polytope([],[],[1,-2],1);
%    cZ2 = cZ1 & P;
%
%    isFullDim(cZ1)
%    isFullDim(cZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper, Adrian Kulmburg
% Written:       02-January-2020 
% Last update:   04-February-2025 (AK, implemented subspace computation)
%                13-February-2025 (TL, added tol)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,2);
tol = setDefaultValues({1e-8},varargin);
inputArgsCheck({ ...
    {cZ,'att','conZonotope'}, ...
    {tol,'att','numeric','scalar'} ...
});

% If the user just wants to know whether the cZ is non-degenerate, this can
% be done quickly

if nargout < 2
    if isempty(cZ.A)
        
        % call zonotope isFullDim method
        res = isFullDim(zonotope(cZ.c,cZ.G),tol);
        
    else
    
        % compute null-space of the constraints
        T = null(cZ.A);
    
        % transform generator matrix into the null-space
        G_ = cZ.G * T;
    
        % check if rank of generator matrix is equal to the dimension
        dimG = size(G_,1);
        rankG = rank(G_,tol);
    
        res = dimG == rankG;
        
    end
else
    % So, the user wants to know the subspace in which the cZ lives.
    % The procedure is very similar to that for zonotopes:
    if isempty(cZ.A)
        
        % call zonotope isFullDim method
        [res, subspace] = isFullDim(zonotope(cZ.c,cZ.G),tol);
        
    else
    
        % compute null-space of the constraints
        T = null(cZ.A);
    
        % transform generator matrix into the null-space
        G_ = cZ.G * T;
    
        % check if rank of generator matrix is equal to the dimension
        dimG = size(G_,1);

        if isempty(G_)
            U = G_;
            Sigma = [];
            r = 0;
        else

            % Since the computation of the rank involves the SVD, and we might
            % need the SVD later on anyhow, we compute it ourselves here:
            [U,Sigma,~] = svd(G_);

            % Very rarely, Sigma can be a vector, then we need to adapt the
            % next line
            if isvector(Sigma)
                s = Sigma(1,1);
            else
                s = diag(Sigma);
            end
            r = sum(s>tol);

        end
    
        res = dimG == r;

        if ~res
            subspace = U(:,1:r);
        else
            subspace = eye(dimG);
        end
        
    end
end

% ------------------------------ END OF CODE ------------------------------
