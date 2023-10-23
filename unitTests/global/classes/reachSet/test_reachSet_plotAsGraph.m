function res = test_reachSet_plotAsGraph
% test_reachSet_plotAsGraph - unit test function for plotting branches of
%    reachSet object
%
% Syntax:
%    res = test_reachSet_plotAsGraph
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
% Written:       07-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% instantiate nx1 reachSet object with different parents
timePoint.set{1} = zonotope(zeros(2,1),eye(2));
timePoint.time{1} = 0;

parent = [0; 1; 1; 2; 3];
R = [];
for i=1:length(parent)
    R = [R; reachSet(timePoint,[],parent(i))];
end

figure;
try

    % no return value
    plotAsGraph(R);
    
    % with return value
    han = plotAsGraph(R);

    % reachable set from continuous-time system (only one branch)
    han = plotAsGraph(R(1));

    close;

catch
    res = false;
    close;
end

% ------------------------------ END OF CODE ------------------------------
