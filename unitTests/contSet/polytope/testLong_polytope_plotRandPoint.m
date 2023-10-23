function res = testLong_polytope_plotRandPoint
% testLong_polytope_plotRandPoint - unit test function of
%    plotRandPoint; this function aims to go through many variations of
%    input arguments
%    note: only run-time errors checked
%
% Syntax:
%    res = testLong_polytope_plotRandPoint
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

% instantiate polytope (via conversion from zonotope)
Z = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
P = polytope(Z);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
%     plotRandPoint(P);
    
    % two arguments: object, dimensions
%     plotRandPoint(P,[1,2]);
%     plotRandPoint(P,[2,3]);

    % three arguments: object, dimensions, number
    plotRandPoint(P,[2,3],50);
    
    % three arguments: object, dimensions, number, linespec
    plotRandPoint(P,[2,3],50,'.b');
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
