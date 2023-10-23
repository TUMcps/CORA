function res = testLong_conPolyZono_plotRandPoint
% testLong_conPolyZono_plotRandPoint - unit test function of
%    plotRandPoint; this function aims to go through many variations of
%    input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLong_conPolyZono_plotRandPoint
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

% construct constrained polynomial zonotope
c = [0;0];
G = [1 0;0 1];
E = [1 0;0 1];
A = [1 -1];
b = 0;
EC = [2 0;0 1];
cPZ = conPolyZono(c,G,E,A,b,EC);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(cPZ);
    
    % two arguments: object, dimensions
    plotRandPoint(cPZ,[1,2]);

    % three arguments: object, dimensions, number
    plotRandPoint(cPZ,[1,2],50);

    % four arguments: object, dimensions, number, linespec
    plotRandPoint(cPZ,[1,2],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
