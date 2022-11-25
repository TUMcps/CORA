function res = testLongDuration_capsule_volume
% testLongDuration_capsule_volume - unit test function of volume
%
% Syntax:  
%    res = testLongDuration_capsule_volume
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
% Written:      28-August-2019
% Last update:  12-March-2021 (MW, add random tests)
% Last revision:---

%------------- BEGIN CODE --------------

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
    tol = 1e-9;
    if abs(vol_true - vol) > tol
        res = false; break;
    end

end


if res
    disp('test_volume successful');
else
    disp('test_volume failed');
end

%------------- END OF CODE --------------