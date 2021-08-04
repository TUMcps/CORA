function res = projectHighDim(obj,N,proj)
% projectHighDim - project a zonotope to a higher-dimensional space
%
% Syntax:  
%    res = projectHighDim(obj,N,dim)
%
% Inputs:
%    obj - zonotope object
%    N - dimension of the higher dimensional space
%    proj - states of the high dimensional space that correspond to the
%          states of the low dimensional mptPolytope object
%
% Outputs:
%    res - zonotope object in the high dimensional space
%
% Example: 
%    Z = zonotope.generateRandom(2);
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

    Z = zeros(N,size(obj.Z,2));
    Z(proj,:) = obj.Z;
    res = zonotope(Z);

%------------- END OF CODE --------------