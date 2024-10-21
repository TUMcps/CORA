function res = isIntersecting_(O,S,type,tol,varargin)
% isIntersecting_ - checks if an empty set intersects with another set
%
% Syntax:
%    res = isIntersecting_(O,S,type,tol)
%
% Inputs:
%    O - emptySet object
%    S - contSet object or numerical vector
%    type - type of check
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    O = emptySet(2);
%    p = [1;1]
%    res = isIntersecting(O,p);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   05-April-2023 (rename isIntersecting_)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% intersection with an empty set is always empty -> always false
res = false;

% ------------------------------ END OF CODE ------------------------------
