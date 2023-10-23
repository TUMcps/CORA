function res = test_capsule_plotRandPoint
% test_capsule_plotRandPoint - unit test function of plotRandPoint; this
%    function aims to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_capsule_plotRandPoint
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

% instantiate capsule
C = capsule([1;1;2],[0;1;-0.5],0.5);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(C);
    
    % two arguments: object, dimensions
    plotRandPoint(C,[1,2]);
    plotRandPoint(C,[2,3]);
    
    % three arguments: object, dimensions, number of points
    plotRandPoint(C,[1,2],50,'.b');
    
    % three arguments: object, dimensions, number of points, linespec
    plotRandPoint(C,[2,3],50,'.b');

    % close figure
    close;
catch ME
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
