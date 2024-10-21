function res = test_matPolytope_printSet
% test_matPolytope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_matPolytope_printSet
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test normal set
V(:,:,1) = [1 2; 0 1];
V(:,:,2) = [1 3; -1 2];
matP = matPolytope(V);

printSet(matP)
printSet(matP,'high')
printSet(matP,'high',true)
printSet(matP,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
