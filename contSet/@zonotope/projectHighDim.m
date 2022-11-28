function Z = projectHighDim(Z,N,proj)
% projectHighDim - project a zonotope to a higher-dimensional space
%
% Syntax:  
%    Z = projectHighDim(Z,N,proj)
%
% Inputs:
%    Z - zonotope object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional mptPolytope object
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
% See also: mptPolytope/projectHighDim

% Author:       Niklas Kochdumper
% Written:      12-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate all-zero center/generator matrix in higher-dimensional space
Znew = zeros(N,size(Z.Z,2));
% insert input argument into desired dimensions of higher-dimensional space
Znew(proj,:) = Z.Z;
% instantiate zonotope object
Z = zonotope(Znew);

%------------- END OF CODE --------------