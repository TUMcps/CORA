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

% test multiple regions ---

% two regions
I1 = interval([1;1],[2;2]);
I2 = interval([-2;-2],[-1;-1]);
pgon = polygon(I1) | polygon(I2);
V_true = [ ...
 NaN, 1, 1, 2, 2, 1, NaN, -2, -2, -1, -1, -2 ; ...
 NaN, 1, 2, 2, 1, 1, NaN, -2, -1, -1, -2, -2 ; ...
];

% plot 
plot(pgon,[1,2])
assert(all(V_true == [ax.Children(1).XData;ax.Children(1).YData] | isnan(V_true),'all'));

% plot with face color (is plotted in multiple regions)
plot(pgon,[1,2],'FaceColor','r')
childs = allchild(ax); % to get children with HandleVisibility off 
assert(compareMatrices(V_true(:,2:6), [childs(2).XData,childs(2).YData]',1e-4,'equal',true));
assert(compareMatrices(V_true(:,8:12), [childs(1).XData,childs(1).YData]',1e-4,'equal',true));

% test holes ---

I1 = interval([-3;-2],[2;3]);
I2 = interval([-1;-1],[1;1]);
pgon = subtract(polygon(I1),polygon(I2));

% plot
plot(pgon);
V = [ax.Children(1).XData;ax.Children(1).YData];
V_true = [ NaN -3 -3 2 2 -3 NaN -1 1 1 -1 -1 ; NaN -2 3 3 -2 -2 NaN -1 -1 1 1 -1 ];
assert(all((V == V_true) | isnan(V_true), 'all'))

% plot
plot(pgon,1:2,'FaceColor','b');
V = [ax.Children(1).XData,ax.Children(1).YData];
V_true = [ -3 -2 ; -3 3 ; 2 3 ; 2 -2 ; 1 -1 ; 1 1 ; -1 1 ; -1 -1 ; 1 -1 ; 2 -2 ];
assert(all((V == V_true), 'all'))

close

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
