function res = test_emptySet_isBounded
% test_emptySet_isBounded - unit test function of isBounded
%
% Syntax:
%    res = test_emptySet_isBounded
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
% See also: -

% Authors:       Tobias Ladner
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% check boundedness
O = emptySet(2);
assert(isBounded(O));

end

% ------------------------------ END OF CODE ------------------------------
