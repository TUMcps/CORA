function res = test_fullspace_polygon
% test_fullspace_polygon - unit test function of polygon
%
% Syntax:
%    res = test_fullspace_polygon
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

% empty case
fs = fullspace.Inf(2);
assertThrowsAs(@polygon, 'CORA:wrongInputInConstructor',fs);

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
