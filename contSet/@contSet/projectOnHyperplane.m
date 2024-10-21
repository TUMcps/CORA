function S_proj = projectOnHyperplane(S,hyp)
% projectOnHyperplane - projects a set onto a hyperplane
%
% Syntax:
%    S_proj = projectOnHyperplane(S,hyp)
%
% Inputs:
%    S - contSet object
%    hyp - polytope object representing a hyperplane
%
% Outputs:
%    S_proj - projected set (of same dimension as input set)
%
% Example: 
%    Z = zonotope([2 1 -1;2 0 1]);
%    hyp = polytope([],[],[1 1],1);
%    res = projectOnHyperplane(Z,hyp);
%
%    figure; hold on; xlim([-3,5]); ylim([-3,5]);
%    plot(hyp,[1,2],'r');
%    plot(Z,[1,2],'g');
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/project

% Authors:       Niklas Kochdumper
% Written:       13-December-2019
% Last update:   24-September-2024 (MW, moved from conHyperplane)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that polytope is a hyperplane
if ~representsa_(hyp,'conHyperplane',1e-12)
    throw(CORAerror('CORA:wrongValue','second','must represent a hyperplane.'));
end

% dimension
n = dim(S);

% normalize hyperplane
hyp_norm = normalizeConstraints(hyp,'A');
c = hyp_norm.Ae'; d = hyp_norm.be;

% linear map A*x + b for the projection
A = eye(n) - c*c';
b = d*c;

% project the set
S_proj = A*S + b;

% ------------------------------ END OF CODE ------------------------------
