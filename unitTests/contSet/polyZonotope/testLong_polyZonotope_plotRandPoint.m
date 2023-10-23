function res = testLong_polyZonotope_plotRandPoint
% testLong_polyZonotope_plotRandPoint - unit test function of
%    plotRandPoint; this function aims to go through many variations of
%    input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLong_polyZonotope_plotRandPoint
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

% instantiate polynomial zonotope
c = rand(4,1)-0.5*ones(4,1);
G = rand(4,6)-0.5*ones(4,6);
ind = datasample(1:6,4,'Replace',false);
G(:,ind) = G(:,ind)./10;
GI = rand(4,2)-0.5*ones(4,2);
E = [eye(4), round(rand(4,2)*5)];
pZ = polyZonotope(c,G,GI,E);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plotRandPoint(pZ);
    
    % two arguments: object, dimensions
    plotRandPoint(pZ,[1,2]);
    plotRandPoint(pZ,[2,3]);
    
    % three arguments: object, dimensions, number
    plotRandPoint(pZ,[1,2],50);
    
    % four arguments: object, dimensions, number, linespec
    plotRandPoint(pZ,[1,2],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
