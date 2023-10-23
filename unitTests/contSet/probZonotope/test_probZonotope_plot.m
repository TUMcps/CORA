function res = test_probZonotope_plot
% test_probZonotope_plot - unit test function of plot; this function aims 
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:
%    res = test_probZonotope_plot
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

% Authors:       Tobias Ladner
% Written:       25-May-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% instantiate probZonotope
Z = [1, 1, -2; 0, 1, 1];
G = [0.6, 1.2; 0.6, -1.2];
S = probZonotope(Z, G);

try
    % try all variations in plotting
    figure; hold on;
    pos = [-57.95, -50.30, 0.76];
    set(gca, 'CameraPosition', pos)
    xlim([-6, 6]);
    ylim([-6, 6]);

    % one argument: object
    plot(S);

    % two arguments: object, dimensions
    plot(S,1);
    plot(S,[1,2]);
    
    % close figure
    close;
catch ME
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
