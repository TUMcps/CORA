function res = test_specification_plotOverTime
% test_specification_plotOverTime - unit test for plotOverTime
%
% Syntax:
%    res = test_specification_plotOverTime
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
% Written:       30-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% init helpers
set1 = zonotope([0;0],[1,-0.7;0.2,1]);
set2 = interval([-1;-2],[2;3]);
time1 = interval(2,4);
time2 = interval(1,3);
spec1 = specification(set1,'unsafeSet',time1);
spec2 = specification(set2,'safeSet',time2);
spec3 = specification(set2,'safeSet');

spec12 = [spec1;spec2];
spec13 = [spec1;spec3];

figure;

try
    % single specification, unsafe
    plotOverTime(spec1);
    hold on

    % single specification, safe
    plotOverTime(spec2,2);
    hold off

    % two specification, both with time
    plotOverTime(spec12);

    % two specifications, only one with time
    plotOverTime(spec13);

    close;

catch ME 
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
