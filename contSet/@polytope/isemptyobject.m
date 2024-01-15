function res = isemptyobject(P)
% isemptyobject - checks if a polytope object is fully empty; this
%    represents R^n since there are no constraints excluding any point
%
% Syntax:
%    res = isemptyobject(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    P = polytope([1 0;-1 0;0 1;0 -1],[3;0;3;-4]);
%    isemptyobject(P); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no inequality or equality constraints
res = isempty(P.b) && isempty(P.be);

% ------------------------------ END OF CODE ------------------------------
