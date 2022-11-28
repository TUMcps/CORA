function res = test_conHyperplane_plot
% test_conHyperplane_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = test_conHyperplane_plot
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

% instantiate constrained hyperplane
c = [1, 1];
d = 2;
A = [1, 0];
b = 2.5;
hyp = conHyperplane(c, d, A, b);

try
    % try all variations in plotting
    figure;
    xlim([-2,2]);
    ylim([-2,2]);
    
    % one argument: object
    plot(hyp);
    
    % two arguments: object, dimensions
    plot(hyp,1);
    plot(hyp,[1,2]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

%------------- END OF CODE --------------
