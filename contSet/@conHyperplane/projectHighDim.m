function res = projectHighDim(hyp,N,dims)
% projectHighDim - projects a constrained hyperplane to a
%    higher-dimensional space
%
% Syntax:  
%    res = projectHighDim(hyp,N,dims)
%
% Inputs:
%    hyp - conHyperplane object
%    N - dimension of the higher-dimensional space
%    dims - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    res - conHyperplane object in the higher-dimensional space
%
% Example: 
%    hs = halfspace([1, -2], 1);
%    C = [-2 -0.5;1 0];
%    d = [-4.25;2.5];
%    hyp = conHyperplane(hs,C,d);
%
%    hyp_ = projectHighDim(hyp,10,[7,9])
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

% check input arguments
inputArgsCheck({{hyp,'att','conHyperplane'};
                {N,'att','numeric',{'nonnan','scalar','nonnegative','integer'}};
                {dims,'att','numeric',{'nonnan','vector','nonnegative'}}});

% project the halfspace
h = projectHighDim(hyp.h,N,dims);

if ~isempty(hyp.C)
    % project the constraints by creating a mptPolytope object
    poly = mptPolytope(hyp.C,hyp.d);
    temp = projectHighDim(poly,N,dims);
    poly = get(temp,'P');

    % construct the resulting high dimensional halfspace object
    res = conHyperplane(h,poly.A,poly.b);
else
    res = conHyperplane(h);
end

%------------- END OF CODE --------------