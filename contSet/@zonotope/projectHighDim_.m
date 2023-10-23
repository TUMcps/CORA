function Z = projectHighDim_(Z,N,proj)
% projectHighDim_ - project a zonotope to a higher-dimensional space
%
% Syntax:
%    Z = projectHighDim_(Z,N,proj)
%
% Inputs:
%    obj - zonotope object
%    N - dimension of the higher dimensional space
%    proj - states of the high dimensional space that correspond to the
%          states of the low dimensional polytope object
%
% Outputs:
%    Z - zonotope object in the higher-dimensional space
%
% Example: 
%    Z = zonotope([-1;1],[3 2 -1; 2 -1 2]);
%    Z_ = projectHighDim(Z,5,[1,3])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/projectHighDim

% Authors:       Niklas Kochdumper
% Written:       12-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate all-zero center/generator matrix in higher-dimensional space
cnew = zeros(N,1);
Gnew = zeros(N,size(Z.G,2));
% insert input argument into desired dimensions of higher-dimensional space
cnew(proj,:) = Z.c;
Gnew(proj,:) = Z.G;
% instantiate zonotope object
Z = zonotope(cnew,Gnew);

% ------------------------------ END OF CODE ------------------------------
