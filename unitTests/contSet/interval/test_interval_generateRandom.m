function res = test_interval_generateRandom
% test_interval_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_interval_generateRandom
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      27-Sep-2019
% Last update:  19-May-2022 (name-value pair syntax)
% Last revision:---

%------------- BEGIN CODE --------------

% empty call
I = interval.generateRandom();

% values for tests
n = 3;

% only dimension
I = interval.generateRandom('Dimension',n);
res = dim(I) == n;


% unify results
res = all(res);

%------------- END OF CODE --------------