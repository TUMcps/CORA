function res = test_fullspace_isBounded
% test_fullspace_isBounded - unit test function of isBounded
%
% Syntax:
%    res = test_fullspace_isBounded
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
fs = fullspace(2);
assert(~isBounded(fs));

% check boundedness of 0-dimensional set
fs = fullspace(0);
assert(isBounded(fs)); % as there is only one element in the set (zeros(0)) 

end

% ------------------------------ END OF CODE ------------------------------
