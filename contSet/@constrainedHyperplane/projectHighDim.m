function res = projectHighDim(obj,N,dim)
% projectHighDim - project a constrainedHyperplane object to a higher
%                  dimensionol space
%
% Syntax:  
%    res = projectHighDim(obj,N,dim)
%
% Inputs:
%    obj - constrainedHyperplane object
%    N - dimension of the higher dimensional space
%    dim - states of the high dimensional space that correspond to the
%          states of the low dimensional mptPolytope object
%
% Outputs:
%    res - constrainedHyperplane object in the high dimensional space
%
% Example: 
%    hs = halfspace([1, -2], 1);
%    C = [-2 -0.5;1 0];
%    d = [-4.25;2.5];
%    ch = constrainedHyperplane(hs,C,d);
%
%    chsHigh = projectHighDim(ch,10,[7,9])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope/projectHighDim, halfspace/projectHighDim

% Author:       Niklas Kochdumper
% Written:      16-July-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % project the halfspace
    h = projectHighDim(obj.h,N,dim);

    % project the constraints by creating a mptPolytope object
    poly = mptPolytope(obj.C,obj.d);
    temp = projectHighDim(poly,N,dim);
    poly = get(temp,'P');

    % construct the resulting high dimensional halfspace object
    res = constrainedHyperplane(h,poly.A,poly.b);

%------------- END OF CODE --------------