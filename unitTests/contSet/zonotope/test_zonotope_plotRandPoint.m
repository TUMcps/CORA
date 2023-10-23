function res = test_zonotope_plotRandPoint
% test_zonotope_plotRandPoint - unit test function of plotRandPoint; this
%    function aims to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_zonotope_plotRandPoint
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

% instantiate zonotope
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(Z);
    
    % two arguments: object, dimensions
    plotRandPoint(Z,[1,2]);
    plotRandPoint(Z,[2,3]);
    
    % three arguments: object, dimensions, number
    plotRandPoint(Z,[2,3],50);
    
    % four arguments: object, dimensions, number, linespec
    plotRandPoint(Z,[2,3],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
