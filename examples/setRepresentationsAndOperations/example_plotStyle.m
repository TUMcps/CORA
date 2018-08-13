function completed = example_plotStyle()
% updated: 21-April-2018, MA

Z1 = zonotope([1 1 1; 1 -1 1]); % create zonotope Z1
Z2 = zonotope([0 1 1; 1 -1 1]); % create zonotope Z1

figure; hold on
plot(Z1); % plot projection of first two coordinate (standard plot)
plot(Z2); % plot projection of first two coordinate (standard plot)
figure; hold on
plot(Z1,[1,2],'r:'); % plot using standard linespec
plot(Z2,[1,2],'r:'); % plot using standard linespec
figure; hold on
plotFilled(Z1,[1,2],'w','EdgeColor','b'); % plot using custom style 'filledFrame'
plotFilled(Z2,[1,2],'w','EdgeColor','b'); % plot using custom style 'filledFrame'
figure; hold on
plotFilled(Z1,[1,2],[.75 .75 .75],'EdgeColor','none'); % plot using custom style 'lightgray'
plotFilled(Z2,[1,2],[.75 .75 .75],'EdgeColor','none'); % plot using custom style 'lightgray'
figure; hold on
plotFilled(Z1,[1,2],[.6 .6 .6],'EdgeColor','none'); % plot using custom style 'darkgray'
plotFilled(Z2,[1,2],[.6 .6 .6],'EdgeColor','none'); % plot using custom style 'darkgray'
figure; hold on
plot(Z1,[1,2],'k-','lineWidth',2); % plot using custom style 'blackEdge'
plot(Z2,[1,2],'k-','lineWidth',2); % plot using custom style 'blackEdge'
figure; hold on
plot(Z1,[1,2],'k-','lineWidth',3); % plot using custom style 'blackEdgeThick'
plot(Z2,[1,2],'k-','lineWidth',3); % plot using custom style 'blackEdgeThick'

%example completed
completed = 1;

%------------- END OF CODE --------------