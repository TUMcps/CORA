function res = testLong_capsule_volume
% testLong_capsule_volume - unit test function of volume
%
% Syntax:
%    res = testLong_capsule_volume
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
% Written:       28-August-2019
% Last update:   12-March-2021 (MW, add random tests)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrOfTests = 1000;
    
for i=1:nrOfTests

    % random dimension
    n = randi([2,10]);

    % instantiate capsule as ball
    c = randn(n,1);
    g = zeros(n,1);
    r = rand(1);
    C = capsule(c,g,r);

    % calculate volume
    vol = volume(C);

    % true volume (n-dim sphere)
    vol_true = (pi^(n/2) / gamma(1+n/2)) * r^n;

    % compare results
    if ~withinTol(vol_true,vol)
        res = false;
        break;
    end

end

% ------------------------------ END OF CODE ------------------------------
