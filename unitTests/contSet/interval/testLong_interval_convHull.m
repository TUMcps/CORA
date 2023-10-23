function res = testLong_interval_convHull
% testLong_interval_convHull - unit test function of convHull
%
% Syntax:
%    res = testLong_interval_convHull
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
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Random cases
res = true;

nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(50);

    % init random interval in [-1,1]
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % shift interval
    I_low = I - rand(n,1);
    I_high = I + rand(n,1);

    % convex hull
    CH = convHull(I_low,I_high);

    % original sets have to be contained
    if ~contains(CH,I_low) || ~contains(CH,I_high)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
