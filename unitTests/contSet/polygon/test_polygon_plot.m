function res = test_polygon_plot()
% test_polygon_plot - unit test function for polygon/plot
%
% Syntax:
%    res = test_polygon_plot()
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
% See also: polygon

% Authors:       Tobias Ladner
% Written:       25-May-2023 (TL, split unit tests)
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate data
x = gallery('uniformdata',30,1,1);
y = gallery('uniformdata',30,1,10);
ind = boundary(x,y);
x = x(ind);
y = y(ind);

% get polygon
pgon = polygon(x,y);

% plot
figure;
plot(pgon)
plot(pgon,[1,2]);

% plot special sets
plot(polygon())
plot(polygon([1 2], [3 4]))

% test colors
plot(pgon,[1,2],'r');
plot(pgon,[1,2],'EdgeColor','r');
plot(pgon,[1,2],'FaceColor','r');
plot(pgon,[1,2],'EdgeColor','k','FaceColor','r');

 % check if plotted correctly
ax = gca();
colorOrder = ax.ColorOrder;

% plot set
plot(pgon,[1,2]);
V = fliplr([x';y']); % somehow get flipped
% check points
assert(compareMatrices(V, [ax.Children(1).XData;ax.Children(1).YData],1e-4,'equal',true));
% test color
assert(isequal(colorOrder(1,:), ax.Children(1).Color));

close

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
