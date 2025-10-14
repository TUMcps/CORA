function res = example_CORAcolor()
% example_CORAcolor - example plotting sets with CORAcolor
%
% Syntax:
%    res = example_CORAcolor()
%
% Inputs:
%    -
%
% Outputs:
%    res - 
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller, Tobias Ladner
% Written:       25-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Create a few sets.
Z = zonotope([0;2],[eye(2) [1;1]]);
P = polytope([-1 3/2 3/2 0; -3 -2 -4 -9/2]);
E = ellipsoid([2 1; 1 2],[-3;1]);

% init figure
figure;

% simple plotting with CORAcolor
subplot(1,2,1); hold on;
title('Simple plotting with CORAcolor')
plot(Z,1:2,'Color',CORAcolor("CORA:blue"))
plot(P,1:2,'Color',CORAcolor("CORA:red"))
plot(E,1:2,'Color',CORAcolor("CORA:yellow"))

% using useCORAcolors
subplot(1,2,2); hold on;
title('Using useCORAcolors')
useCORAcolors("CORA:default"); % defines the color order
plot(Z)
plot(P)
plot(E)

disp('See also <a href="matlab:open CORA_PLOT_FILLED">CORA_PLOT_FILLED</a> to toggle plotting sets with face color per default.')

% example completed
res = true;

end


% ------------------------------ END OF CODE ------------------------------
