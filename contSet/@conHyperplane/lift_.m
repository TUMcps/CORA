function res = lift_(hyp,N,dims)
% lift_ - lifts a constrained hyperplane to a higher-dimensional space
%
% Syntax:
%    res = lift_(hyp,N,dims)
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
%    hyp_ = lift(hyp,10,[7,9])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, polytope/lift, halfspace/lift

% Authors:       Niklas Kochdumper
% Written:       16-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% project the halfspace
h = lift_(hyp.h,N,dims);

if ~isempty(hyp.C)
    % project the constraints by creating a polytope object
    P = polytope(hyp.C,hyp.d);
    P = lift_(P,N,dims);
    
    % construct the resulting high dimensional halfspace object
    res = conHyperplane(h,P.A,P.b);
else
    res = conHyperplane(h);
end

% ------------------------------ END OF CODE ------------------------------
