function res = test_fullspace_interval
% test_fullspace_interval - unit test function of interval
%
% Syntax:
%    res = test_fullspace_interval
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

% convert to interval
I = interval(fs);

% true result
I_true = interval(-Inf(n,1),Inf(n,1));

% compare results
res = I == I_true;

% ------------------------------ END OF CODE ------------------------------
