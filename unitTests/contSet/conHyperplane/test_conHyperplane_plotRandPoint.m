function res = test_conHyperplane_plotRandPoint
% test_conHyperplane_plotRandPoint - unit test function of plotRandPoint;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_conHyperplane_plotRandPoint
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
% Written:       11-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% instantiate constrained hyperplane
a = [1, 1, 1]; b = 2;
hyp = conHyperplane(a, b);

try
    % try all variations in plotting
    figure;
    xlim([-2,2]);
    ylim([-2,2]);
    
    % one argument: object
    plotRandPoint(hyp);
    
    % two arguments: object, dimensions
    plotRandPoint(hyp,[1,2]);
    plotRandPoint(hyp,[2,3]);

    % three arguments: object, dimensions, number
    plotRandPoint(hyp,[1,2],50);

    % three arguments: object, dimensions, number
    plotRandPoint(hyp,[1,2],75,'.b');
    
    % close figure
    close;
catch ME
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
