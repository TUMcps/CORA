function res = test_conPolyZono_plot
% test_conPolyZono_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_conPolyZono_plot
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

% Author:       Mark Wetzlinger
% Written:      25-May-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% construct constrained polynomial zonotope
c = [0;0];
G = [1 0;0 1];
expMat = [1 0;0 1];
A = [1 -1];
b = 0;
expMat_ = [2 0;0 1];
cPZ = conPolyZono(c,G,expMat,A,b,expMat_);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(cPZ);
    
    % two arguments: object, dimensions
    plot(cPZ,1);
    plot(cPZ,[1,2]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------