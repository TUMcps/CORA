function res = testLong_interval_volume
% testLong_interval_volume - unit test function of volume computation
%
% Syntax:
%    res = testLong_interval_volume
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
    lowerLimits = -3+3*rand(n,1);
    upperLimits = 3*rand(n,1);
    I = interval(lowerLimits, upperLimits);
    
    % compute volume
    vol = volume(I);
    
    % true volume
    vol_true = prod(upperLimits - lowerLimits);
    
    % compare results
    if ~withinTol(vol,vol_true,tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
