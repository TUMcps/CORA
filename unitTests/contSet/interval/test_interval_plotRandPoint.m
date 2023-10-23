function res = test_interval_plotRandPoint
% test_interval_plotRandPoint - unit test function of plotRandPoint; this
%    function aims to go through many variations of input arguments
%    note: only run-time errors checked
%
% Syntax:
%    res = test_interval_plotRandPoint
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

% init interval
I = interval([1;1;2],[3;4;7]);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(I);
    
    % two arguments: object, dimensions
    plotRandPoint(I,[1,2]);
    plotRandPoint(I,[2,3]);
    
    % three arguments: object, dimensions, number
    plotRandPoint(I,[1,2],50);
    
    % three arguments: object, dimensions, number, linespec
    plotRandPoint(I,[1,2],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
