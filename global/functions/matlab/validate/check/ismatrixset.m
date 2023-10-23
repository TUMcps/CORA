function res = ismatrixset(S)
% ismatrixset - checks whether a set is a matrix set (object of either
%    intervalMatrix, matPolytope, or matZonotope classes); mainly used to
%    avoid clutter in case differentations
%
% Syntax:
%    res = ismatrixset(S)
%
% Inputs:
%    S - anything
%
% Outputs:
%    res - true/false
%
% Example:
%    C = capsule([1;1],[1;-1],0.5);
%    matP = matPolytope({[1 2; 0 1],[1 3; -1 2]});
%
%    ismatrixset(C)
%    ismatrixset(matP)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       17-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;

if isa(S,'intervalMatrix') || isa(S,'matPolytope') || isa(S,'matZonotope')
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
