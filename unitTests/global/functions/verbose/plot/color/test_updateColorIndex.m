function res = test_updateColorIndex
% test_updateColorIndex - unit test function for updating the color index
%
% Syntax:  
%    res = test_updateColorIndex()
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
% See also: colororder

% Author:       Tobias Ladner
% Written:      13-July-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = [];

% create figure
f = figure;
ax = gca();

% test with hold off
hold off
oldColorIndex = 1;
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;

plotNPlots(1);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;

plotNPlots(5);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;


% test with hold on
hold on
oldColorIndex = ax.ColorOrderIndex;
plotNPlots(1);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

oldColorIndex = ax.ColorOrderIndex;
plotNPlots(2);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

oldColorIndex = ax.ColorOrderIndex;
plotNPlots(5);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

close(f);

% gather results
res = all(resvec);

end

% Auxiliary functions -----------------------------------------------------

function plotNPlots(n)
    color = CORAcolor('CORA:next');
    for i=1:n
        plot(1:5, rand(5, 1),'Color',color);
    end
end


%------------- END OF CODE --------------
