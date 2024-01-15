function S = projectOnHyperplane(hyp,S)
% projectOnHyperplane - projects a set onto a hyperplane
%
% Syntax:
%    S = projectOnHyperplane(hyp, S)
%
% Inputs:
%    hyp - conHyperplane object
%    S - contSet object
%
% Outputs:
%    S - projected set
%
% Example: 
%    hyp = conHyperplane([1 1],1);
%    Z = zonotope([2 1 -1;2 0 1]);
%    res = projectOnHyperplane(hyp,Z);
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get object properties
c = hyp.a'; d = hyp.b;

% normalization
temp = norm(c);

c = c./temp;
d = d/temp;

% linear map A*x + b for the projection
A = eye(length(c)) - c*c';
b = d*c;

% project the set
S = A*S + b;

% ------------------------------ END OF CODE ------------------------------
