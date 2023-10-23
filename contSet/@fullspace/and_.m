function S = and_(fs,S,varargin)
% and_ - overloads '&' operator, computes the intersection of a
%    full-dimensional space and another set or numerical vector;
%    case R^0: can only intersect with R^0, 0 (not representable in
%    MATLAB), or the empty set, resulting in R^0, R^0, and the empty set,
%    respectively
%
% Syntax:
%    S = and_(fs,S)
%    S = and_(fs,S,method)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object or numerical vector
%    method - (optional) approximation method
%
% Outputs:
%    set - intersection
%
% Example: 
%    fs = fullspace(2);
%    S = zonotope([1;1],[2 1; -3 1]);
%    S = fs & S;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (rename and_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% re-order sets due to class preference
[fs,S] = findClassArg(fs,S,'fullspace');

% intersection is always the other set

% ------------------------------ END OF CODE ------------------------------
