function res = test_halfspace_plot
% test_halfspace_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_halfspace_plot
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

% normal vector and offset
a = [-0.5; -1; 0.1];
b = -1;

% instantiate constrained hyperplane
hs = halfspace(a,b);

try
    % try all variations in plotting
    figure;
    xlim([-2,2]);
    ylim([-2,2]);
    
    % one argument: object
    plot(hs);
    
    % two arguments: object, dimensions
    plot(hs,1);
    plot(hs,[1,2]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------
