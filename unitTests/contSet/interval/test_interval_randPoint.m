function res = test_interval_randPoint
% test_interval_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_interval_randPoint
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

% Author:       Mark Wetzlinger
% Written:      21-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
I = interval();
res = isempty(randPoint(I));

%------------- END OF CODE --------------