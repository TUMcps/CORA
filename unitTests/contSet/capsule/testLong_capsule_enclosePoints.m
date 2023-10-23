function res = testLong_capsule_enclosePoints
% testLong_capsule_enclosePoints - unit test function of enclosePoints
%
% Syntax:
%    res = testLong_capsule_enclosePoints
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrTests = 1000;

for i=1:nrTests
    % random dimension
    n = randi([2,8]);

    % random number of points
    nrPoints = randi(100);

    % init points
    p = randn(n,nrPoints);

    % compute enclosing capsule
    C = capsule.enclosePoints(p);

    % check if all points are contained in capsule
    if ~all(contains(C,p))
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
