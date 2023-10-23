function completed = example_plot_posXYZ()
% example_plot_posXYZ - shows how sets can be plotted in the XYZ space
%
% Syntax:
%    completed = example_plot_posXYZ()
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

figure

% plotting high-dimensional intervals -------------------------------------

subplot(2, 2, 1); hold on;

I = interval([3;4;2;7;4;5],[4;4;7;8;5;7]);
for i=1:dim(I)
    plot(I, i, 'YPos', i); % plots dimension i to y=i
end

subplot(2, 2, 2); hold on;
for i=1:dim(I)
    plot(I, i, 'XPos', i); % plots dimension i to x=i
end

% plotting in z-direction -------------------------------------------------

subplot(2, 1, 2); hold on;
Z0 = zonotope([1;1], 0.2*[0 1 1; 1 1, 0]);
a = pi/10;
M = [cos(a) sin(a);-sin(a) cos(a)];
view(140, 40)

Z = Z0;
for i=1:21
    plot(Z, [1,2], 'ZPos', i, 'FaceColor',colorblind('b'),'EdgeColor','k');
    Z = M*Z;
end

Z = -1 * Z0;
for i=1:21
    plot(Z, [1,2], 'ZPos', i, 'FaceColor',colorblind('r'),'EdgeColor','k');
    Z = M*Z;
end

completed = true;

% ------------------------------ END OF CODE ------------------------------
