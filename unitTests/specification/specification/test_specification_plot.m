function res = test_specification_plot
% test_specification_plot - unit test for plot
%
% Syntax:
%    res = test_specification_plot
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
set = zonotope([0;0],[1,-0.7;0.2,1]);
set2 = interval([-1;-2],[2;3]);
list = {set,set2};
time = interval(2,4);
loc = [1 2];

figure;

try
    % unsafe set
    spec = specification(set,'unsafeSet');
    plot(spec);
    hold on

    % safe set
    spec = specification(set2,'safeSet');
    plot(spec);
    hold off

    % list of sets
    spec = specification(list,'unsafeSet');
    plot(spec);

    % invariant
    spec = specification(set2,'invariant');
    plot(spec);

    % mixed
    spec = [specification(set,'safeSet'),specification(set,'safeSet'),specification(set2,'unsafeSet')];
    plot(spec)

    % with time
    spec = specification(set,'safeSet',time);
    plot(spec);

    % with time and location
    spec = specification(set,'safeSet',time,loc);
    plot(spec);

    close;

catch ME
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
