function res = testLong_conZonotope_plotRandPoint
% testLong_conZonotope_plotRandPoint - unit test function of
%    plotRandPoint; this function aims to go through many variations of
%    input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLong_conZonotope_plotRandPoint
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

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZ = conZonotope(Z,A,b);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(cZ);
    
    % two arguments: object, dimensions
    plotRandPoint(cZ,[1,2]);
    
    % three arguments: object, dimensions, number
    plotRandPoint(cZ,[1,2],50);
    
    % four arguments: object, dimensions, number, linespec
    plotRandPoint(cZ,[1,2],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
