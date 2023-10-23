function res = contains_(hyp,S,varargin)
% contains_ - determines if a constrained hyperplane contains a set or a
%    point
%
% Syntax:
%    res = contains_(hyp,S)
%
% Inputs:
%    hyp - conHyperplane object
%    S - contSet object or single point
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    point = [0;0];
%    contains(hyp,point)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, conHyperplane/isIntersecting_

% Authors:       Victor Gassmann
% Written:       19-July-2021
% Last update:   25-November-2022 (MW, rename 'contains')
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% containment check (input arguments are checked there)
res = contains_(polytope(hyp),S,varargin{:});

% ------------------------------ END OF CODE ------------------------------
