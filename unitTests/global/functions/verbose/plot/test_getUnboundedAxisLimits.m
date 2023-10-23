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

resvec = [];

try
    figure;

    % empty plot
    [xLim,yLim] = getUnboundedAxisLimits();
    resvec(end+1) = all([0 1 0 1] == [xLim yLim]);

    % give weird vertices
    [xLim,yLim] = getUnboundedAxisLimits([1.25, 1.223; 1.34, 1.2]);
    resvec(end+1) = all(withinTol([0 1.4 0 1.4], [xLim yLim]));

    % set axis
    xlim([-2,1]); ylim([1,3]);

    [xLim,yLim] = getUnboundedAxisLimits();
    resvec(end+1) = all(withinTol([-2 1 1 3], [xLim yLim]));

    % with vertices
    [xLim,yLim] = getUnboundedAxisLimits([1.25, 1.223; 1.34, 1.2]);
    resvec(end+1) = all(withinTol([-2 1.5 1 3], [xLim yLim]));

    % zlim ---
    view(3); 
    xlim([-2,1]); ylim([1,3]); zlim([4,5]);
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    resvec(end+1) = all(withinTol([-2 1 1 3 4 5], [xLim yLim zLim]));

    [xLim,yLim,zLim] = getUnboundedAxisLimits([1.25 1.223; 1.34 1.2; 2.2 3]);
    resvec(end+1) = all(withinTol([-2 2 1 3 2 5], [xLim yLim zLim]));

    close;

catch ME
    close
    resvec(end+1) = false;
end

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
