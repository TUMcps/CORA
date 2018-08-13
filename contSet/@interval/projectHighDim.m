function res = projectHighDim(obj,N,dim)
% projectHighDim - project an interval object to a higher dimensional
%                  space
%
% Syntax:  
%    res = projectHighDim(obj,N,dim)
%
% Inputs:
%    obj - interval object
%    N - dimension of the higher dimensional space
%    dim - states of the high dimensional space that correspond to the
%          states of the low dimensional interval object
%
% Outputs:
%    res - interval object in the high dimensional space
%
% Example: 
%    inter = interval([1;-3],[4;1]);
%    interHigh = projectHighDim(inter,10,[3;5]);
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
    sup = ones(N,1) * inf;
    infi = ones(N,1) * (-inf);
    
    % insert parameters from the low dimensional interval object
    sup(dim) = supremum(obj);
    infi(dim) = infimum(obj);

    % construct the resulting high dimensional interval object
    res = interval(infi,sup);

%------------- END OF CODE --------------