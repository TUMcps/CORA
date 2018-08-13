function res = projectHighDim(obj,N,dim)
% projectHighDim - project a halfspace object to a higher dimensional
%                  space
%
% Syntax:  
%    res = projectHighDim(obj,N,dim)
%
% Inputs:
%    obj - halfspace object
%    N - dimension of the higher dimensional space
%    dim - states of the high dimensional space that correspond to the
%          states of the low dimensional mptPolytope object
%
% Outputs:
%    res - halfspace object in the high dimensional space
%
% Example: 
%    C = [2.1, 3.4, 5.2];
%    d = 1.7;
%    hs = halfspace(C,d);
%
%    hsHigh = projectHighDim(hs,10,[1,7,9])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope/projectHighDim, interval/projectHighDim

% Author:       Niklas Kochdumper
% Written:      16-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % initialize variables
    C = zeros(N,1);
    d = obj.d;

    % insert parameters from the original halfspace object
    C(dim) = obj.c;

    % construct the resulting high dimensional halfspace object
    res = halfspace(C,d);

%------------- END OF CODE --------------