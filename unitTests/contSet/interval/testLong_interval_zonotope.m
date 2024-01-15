function res = testLong_interval_zonotope
% testLong_interval_zonotope - unit test function of zonotope conversion
%
% Syntax:
%    res = testLong_interval_zonotope
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
% Written:       04-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
tol = 1e-9;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    
    % create random interval
    n = randi(10);
    lb = -3+3*rand(n,1);
    ub = 3*rand(n,1);
    I = interval(lb, ub);
    
    % convert to zonotope
    Z = zonotope(I);
    
    % true zonotope
    c_true = center(I);
    G_true = diag(rad(I));
    
    % compare results
    if ~all(withinTol(center(Z),c_true,tol))
        throw(CORAerror('CORA:testFailed'));
    elseif ~compareMatrices(generators(Z),G_true,tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
