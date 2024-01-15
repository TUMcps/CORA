function res = testLong_interval_enlarge
% testLong_interval_enlarge - unit test function of enlarge
%
% Syntax:
%    res = testLong_interval_enlarge
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

% Authors:       Mark Wetzlinger
% Written:       03-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% tolerance
tol = 1e-10;

% number of tests
nrTests = 100;

for i=1:nrTests
    % random dimension
    n = floor(1 + 9*rand(1));

    % random bounds
    lb = -10*rand(n,1);
    ub = 10*rand(n,1);
    I = interval(lb, ub);

    % random factor
    factor = 2;

    % enlarge
    I_enlarge = enlarge(I, factor);

    % true solution
    c = (ub + lb)/2;
    r = (ub - lb)/2;
    lower_true = c - r*factor;
    upper_true = c + r*factor;
    I_true = interval(lower_true, upper_true);

    % check
    if ~isequal(I_enlarge,I_true,tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
