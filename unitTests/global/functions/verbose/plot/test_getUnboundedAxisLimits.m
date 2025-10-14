function res = test_getUnboundedAxisLimits
% test_getUnboundedAxisLimits - unit test function for getUnboundedAxisLimits
%
% Syntax:
%    res = test_getUnboundedAxisLimits
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

% Authors:       Tobias Ladner
% Written:       26-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

try
    figure;

    % empty plot
    [xLim,yLim] = getUnboundedAxisLimits();
    assert(all([0 1 0 1] == [xLim yLim]));

    % give weird vertices
    V = [1.25, 1.223; 1.34, 1.2];
    [xLim,yLim] = getUnboundedAxisLimits(V);
    assert(all([xLim(1);yLim(1)] <= V | V <= [xLim(2);yLim(2)],"all"));

    % set axis
    xlim([-2,1]); ylim([1,3]);

    [xLim,yLim] = getUnboundedAxisLimits();
    assert(all(withinTol([-2 1 1 3], [xLim yLim])));

    % with vertices
    V = [1.25, 1.223; 1.34, 1.2];
    [xLim,yLim] = getUnboundedAxisLimits(V);
    assert(all([xLim(1);yLim(1)] <= V | V <= [xLim(2);yLim(2)],"all"));

    % zlim ---
    view(3); 
    xlim([-2,1]); ylim([1,3]); zlim([4,5]);
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    assert(all(withinTol([-2 1 1 3 4 5], [xLim yLim zLim])));

    % with vertices
    V = [1.25 1.223; 1.34 1.2; 2.2 3];
    [xLim,yLim,zLim] = getUnboundedAxisLimits(V);
    assert(all([xLim(1);yLim(1);zLim(1)] <= V | V <= [xLim(2);yLim(2);zLim(2)],"all"));

    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
