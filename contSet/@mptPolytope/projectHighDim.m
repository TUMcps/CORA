function res = projectHighDim(obj,N,dim)
% projectHighDim - project a mptPolytope object to a higher dimensional
%                  space
%
% Syntax:  
%    res = projectHighDim(obj,N,dim)
%
% Inputs:
%    obj - mptPolytope object
%    N - dimension of the higher dimensional space
%    dim - states of the high dimensional space that correspond to the
%          states of the low dimensional mptPolytope object
%
% Outputs:
%    res - mptPolytope object in the high dimensional space
%
% Example: 
%    A = [1 0;-1 0;0 1;0 -1;1 1];
%    b = [1;1;1;1;1];
%    poly = mptPolytope(A,b);
%
%    polyHigh = projectHighDim(poly,10,[4,5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      06-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % initialize variables
    A = zeros(size(obj.P.A,1),N);
    b = obj.P.b;

    % insert parameters from the original mptPolytope object
    A(:,dim) = obj.P.A;

    % construct the resulting high dimensional mptPolytope object
    res = mptPolytope(A,b);

%------------- END OF CODE --------------