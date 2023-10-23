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

% Authors:       Tobias Ladner
% Written:       28-February-2023
% Last update:   24-March-2023
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% specify color via LineSpec
% see also https://www.mathworks.com/help/matlab/ref/plot.html#btzitot_sep_mw_3a76f056-2882-44d7-8e73-c695c0c54ca8
plot(Z, [1 2], 'g--', 'DisplayName', 'Green')
Z = M * Z;

% or using EdgeColor
plot(Z, [1 2], 'EdgeColor', [0 0 1], 'DisplayName', 'EdgeColor')
Z = M * Z;

% or fill using FaceColor
plot(Z, [1 2], 'FaceColor', [1 0 0], 'DisplayName', 'FaceColor')
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
for i=1:3
    plot(Z, [1 2], 'DisplayName', 'Default')
    Z = M * Z;
end

% or explicitly specify the EdgeColor according to colororder
for i=1:3
    plot(Z, [1 2], 'EdgeColor', 'next', ...
        'DisplayName', 'EdgeColor')
    Z = M * Z;
end


% or fill with 'FaceColor' and 'next'
for i=1:3
    plot(Z, [1 2], 'FaceColor', 'next', ...
        'DisplayName', 'FaceColor')
    Z = M * Z;
end

% back to default
for i=1:5
    plot(Z, [1 2], 'DisplayName', 'Default')
    Z = M * Z;
end

legend('Location', 'eastoutside')


% ------------------------------ END OF CODE ------------------------------
