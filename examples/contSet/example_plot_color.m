function example_plot_color()
% example_plot_color - shows how the color of the plots can be configured
%
% Syntax:  
%    completed = example_plot_color()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: colororder

% Author:        Tobias Ladner
% Written:       28-February-2023
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% init
Z0 = zonotope([-3;3], [1 0 1; 0 1 1]/5);
a = 0.8;
M = [cos(a) sin(a); -sin(a) cos(a)] * 0.8;

% -------------------------------------------------------------------------
% plotting options
% see also https://de.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html

figure;

subplot(2, 1, 1); hold on;
% default:
Z = Z0;
for i=1:4
    plot(Z, [1 2], 'DisplayName', 'Default')
    Z = M * Z;
end

% specify color
plot(Z, [1 2], 'g', 'DisplayName', 'Green')
Z = M * Z;

% or fill
plot(Z, [1 2], 'EdgeColor', 'k', 'FaceColor', [0 0.4470 0.7410], ...
    'DisplayName', 'Filled')
Z = M * Z;

for i=1:7
    plot(Z, [1 2], 'DisplayName', 'Default')
    Z = M * Z;
end
legend('Location', 'eastoutside')

% -------------------------------------------------------------------------
% change default plotting colors
% see also https://de.mathworks.com/help/matlab/ref/colororder.html

subplot(2, 1, 2); hold on;

% change default plotting colors for current axis
newcolors = [0 0.5 1; 0.5 0 1; 0.7 0.7 0.7];
colororder(newcolors)

Z = Z0;
for i=1:13
    plot(Z, [1 2], 'DisplayName', 'Default')
    Z = M * Z;
end

legend('Location', 'eastoutside')


%------------- END OF CODE --------------
