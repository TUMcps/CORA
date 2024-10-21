function res = test_specification_printSpec
% test_specification_printSpec - unit test function of printSpec
%
% Syntax:
%    res = test_specification_printSpec
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

% test print of a simple system
set = polytope([1 1]/sqrt(2),1);
spec = specification(set,'unsafeSet',interval(2,3));

printSpec(spec)
printSpec(spec,'high')
printSpec(spec,'high',true)
printSpec(spec,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
