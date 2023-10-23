function example_manual_polyZonotope_construction()
% example_manual_polyZonotope_construction - example from the manual 
% demonstrating the construction of a polyZonotope
%
% Syntax:
%   example_manual_polyZonotope_construction()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       29-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct polynomial zonotope
c = [4;4];
G = [2 1 2; 0 2 2];
E = [1 0 3;0 1 1];
GI = [1;0];

pZ = polyZonotope(c,G,GI,E);

% plot --------------------------------------------------------------------

figure;
drawArrow = @(pos, dir) quiver(pos(1),pos(2),dir(1),dir(2), 'k');

% (d): first plot end result to get axis limits
subplot(1,4,4); hold on; box on;
plot(pZ)
scatter(c(1,:),c(2,:),'.k')
enlargeAxis(1.2);
title('(d)')

% read axis limits
xLim = xlim(); yLim = ylim();

% (a)
subplot(1,4,1); hold on; box on;
scatter(c(1,:),c(2,:),'.k')
drawArrow(c,G(:,1)); drawArrow(c,-G(:,1));
drawArrow(c,G(:,2)); drawArrow(c,-G(:,2));
title('(a)')
xlim(xLim); ylim(yLim);

% (b)
subplot(1,4,2); hold on; box on;
plot(polyZonotope(c,G(:,1:2),[],E(:,1:2)))
scatter(c(1,:),c(2,:),'.k')
drawArrow([7;6],G(:,3)); drawArrow([1;2],G(:,3));
drawArrow([3;6],-G(:,3)); drawArrow([5;2],-G(:,3));
title('(b)')
xlim(xLim); ylim(yLim);

% (c)
subplot(1,4,3); hold on; box on;
plot(polyZonotope(c,G,[],E))
scatter(c(1,:),c(2,:),'.k')
drawArrow([1;4],GI(:,1)); drawArrow([1;4],-GI(:,1));
drawArrow([9;8],GI(:,1)); drawArrow([9;8],-GI(:,1));
drawArrow([3.01697;0.0169735],GI(:,1)); drawArrow([3.01697;0.0169735],-GI(:,1));
title('(c)')
xlim(xLim); ylim(yLim);

% ------------------------------ END OF CODE ------------------------------
