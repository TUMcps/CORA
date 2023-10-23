function res = test_ellipsoid_plotRandPoint
% test_ellipsoid_plotRandPoint - unit test function of plotRandPoint
%
% Syntax:
%    res = test_ellipsoid_plotRandPoint
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
load cases.mat E_c

for i=1:length(E_c)
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.Ed1;
    E0 = E_c{i}.E0;
    
    res = aux_tryPlot(E1) && aux_tryPlot(Ed1) && aux_tryPlot(E0);
    
    if ~res
        break;
    end 
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_tryPlot(E)

res = true;
try
    % try all variations in plotting
    h = figure;
    
    % one argument: object
    plotRandPoint(E);

    % two arguments: object and dimension
    plotRandPoint(E,[1,2]);
    
    % three arguments: object, dimensions, number
    plotRandPoint(E,[1,2],50);
    
    % four arguments: object, dimensions, number, linespec, NVpairs
    plotRandPoint(E,[1,2],50,'.b');
    
    % close figure
    close(h);
catch
    close;
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
