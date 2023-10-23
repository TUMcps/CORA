function res = isFullDim(cZ)
% isFullDim - checks if the dimension of the affine hull of a constrained
%    zonotope is equal to the dimension of its ambient space
%
% Syntax:
%    res = isFullDim(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%
%    hyp = conHyperplane([1,-2],1);
%    cZ2 = cZ1 & hyp;
%
%    isFullDim(cZ1)
%    isFullDim(cZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Authors:       Niklas Kochdumper
% Written:       02-January-2020 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isempty(cZ.A)
    
    % call zonotope isFullDim method
    res = isFullDim(zonotope(cZ.c,cZ.G));
    
else

    % compute null-space of the constraints
    T = null(cZ.A);

    % transform generator matrix into the null-space
    G_ = cZ.G * T;

    % check if rank of generator matrix is equal to the dimension
    dimG = size(G_,1);
    rankG = rank(G_);

    res = dimG == rankG;
    
end

% ------------------------------ END OF CODE ------------------------------
