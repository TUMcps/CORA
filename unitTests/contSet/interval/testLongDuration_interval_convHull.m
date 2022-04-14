function res = testLongDuration_interval_convHull
% testLongDuration_interval_convHull - unit test function of convHull
%
% Syntax:
%    res = testLongDuration_interval_convHull
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
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Random cases
res_rand = true;

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
    if ~in(CH,I_low) || ~in(CH,I_high)
        res_rand = false; break;
    end

end

% dimension mismatch
res_mismatch = true;
I_dim3 = interval(-rand(3,1),rand(3,1));
I_dim7 = interval(-rand(7,1),rand(7,1));
try
    convHull(I_dim3,I_dim7);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res_mismatch = false;
    end
end

% combine results
res = res_rand;

if res
    disp('test_convHull successful');
else
    disp('test_convHull failed');
end

%------------- END OF CODE --------------

