function res = test_fullspace_randPoint
% test_fullspace_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_fullspace_randPoint
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

% Authors:       Mark Wetzlinger
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 3;
fs = fullspace(n);

% sample point
p = randPoint(fs);
res = contains(fs,p);

% more syntax variations
p = randPoint(fs,5);
res(end+1,1) = all(contains(fs,p));
p = randPoint(fs,10,'standard');
res(end+1,1) = all(contains(fs,p));

% extreme points
p = randPoint(fs,1,'extreme');
res(end+1,1) = all(contains(fs,p));
p = randPoint(fs,'all','extreme');
res(end+1,1) = all(contains(fs,p));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
