function res = test_zonoBundle_plotRandPoint
% test_zonoBundle_plotRandPoint - unit test function of plotRandPoint; this
%    function aims to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_zonoBundle_plotRandPoint
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

% instantiate zonotope bundle
Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB = zonoBundle({Z1,Z2});

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(zB);

    % two arguments: object, dimensions
    plotRandPoint(zB,[1,2]);
    plotRandPoint(zB,[2,3]);

    % three arguments: object, dimensions, number
    plotRandPoint(zB,[2,3],50);

    % four arguments: object, dimensions, number, linespec
    plotRandPoint(zB,[2,3],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
