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

% Authors:       Tobias Ladner
% Written:       13-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% create figure
f = figure;
ax = gca();

% test with hold off
hold off
oldColorIndex = 1;
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;

aux_plotNPlots(1);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;

aux_plotNPlots(5);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == 2;


% test with hold on
hold on
oldColorIndex = ax.ColorOrderIndex;
aux_plotNPlots(1);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

oldColorIndex = ax.ColorOrderIndex;
aux_plotNPlots(2);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

oldColorIndex = ax.ColorOrderIndex;
aux_plotNPlots(5);
updateColorIndex(oldColorIndex)
resvec(end+1) = ax.ColorOrderIndex == (oldColorIndex+1);

% reset color index (pass index < 1)
updateColorIndex(0)
resvec(end+1) = ax.ColorOrderIndex == 1;

updateColorIndex(-1)
resvec(end+1) = ax.ColorOrderIndex == 1;

close(f);

% gather results
res = all(resvec);

end


% Auxiliary functions -----------------------------------------------------

function aux_plotNPlots(n)
    color = CORAcolor('CORA:next');
    for i=1:n
        plot(1:5, rand(5, 1),'Color',color);
    end
end


% ------------------------------ END OF CODE ------------------------------
