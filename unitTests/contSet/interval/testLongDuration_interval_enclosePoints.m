function res = testLongDuration_interval_enclosePoints
% testLongDuration_interval_enclosePoints - unit test function of
%    enclosePoints
%
% Syntax:  
%    res = testLongDuration_interval_enclosePoints
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

% random cases
res = true;
nrOfTests = 1000;
nrPts = 100;

for i=1:nrOfTests

    % random dimension
    n = randi(50);
    
    % bounds for points
    minbound = -5*rand(n,1);
    maxbound = 5*rand(n,1);
    % diff
    diameter = maxbound - minbound;

    % generate random point clouds (from -1 to 1)
    pts = minbound + diameter .* rand(n,nrPts);

    % init enclosing interval (has to be smaller than -1 to 1)
    I = interval.enclosePoints(pts);

    % check with correct solution
    if any(infimum(I) < minbound) || any(supremum(I) > maxbound)
        res = false; break;
    end

end


if res
    disp('test_enclosePoints successful');
else
    disp('test_enclosePoints failed');
end

%------------- END OF CODE --------------
