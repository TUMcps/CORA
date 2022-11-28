function res = test_zonoBundle_plot
% test_zonoBundle_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:  
%    res = test_zonoBundle_plot
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
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

% instantiate zonotope bundle
Z1 = zonotope([1;-1;2],[2 -1 3; 0 1 -1; -1 4 2]);
Z2 = Z1 + [1;0;0];
zB = zonoBundle({Z1,Z2});

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(zB);

    % two arguments: object, dimensions
    plot(zB,1);
    plot(zB,[2,3]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------