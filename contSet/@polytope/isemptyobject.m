function res = isemptyobject(P)
% isemptyobject - checks if a polytope object is fully empty;
%    if the H representation is given, this represents R^n since there are
%    no constraints excluding any point
%    if, instead, the V representation is given, this represents the empty
%    set since there are no vertices
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
% Last update:   14-July-2024 (MW, support V representation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no inequality or equality constraints
res_H = P.isHRep.val && isempty(P.b_.val) && isempty(P.be_.val);
% no vertices
res_V = P.isVRep.val && isempty(P.V_.val);

% combine information
res = res_H || res_V;

% ------------------------------ END OF CODE ------------------------------
