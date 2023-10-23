function res = testLong_levelSet_plot
% testLong_levelSet_plot - unit test function of plot
%
% Syntax:
%    res = testLong_levelSet_plot
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
% See also: plot.m

% Authors:       Maximilian Perschl
% Written:       08-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% define level sets
syms x y z

% 2d equality:
eq = x^2 + y^2 - 4;
levelSets{1} = levelSet(eq,[x;y],'==');
% 2d inequality:
eq = x^2 + y^2 -4;
levelSets{2} = levelSet(eq,[x;y],'<=');
% 3d equality non-solvable:
eq = x^2 + y^2 + z^2 - 4;
levelSets{3} = levelSet(eq,[x;y;z],'==');
% 3d equality solvable:
eq = x + y^2 + z^2 - 4;
levelSets{4} = levelSet(eq,[x;y;z],'==');
% 3d inequality:
eq = x^2 + y^2 + z^2 - 4;
levelSets{5} = levelSet(eq,[x;y;z],'<=');
% 1d plot:
eq = x^2 - 4;
levelSets{6} = levelSet(eq,x,'==');

% loop over sets
for i = 1:length(levelSets)
    figure;
    xlim([-3,3]);
    ylim([-3,3]);

    try
        if levelSets{i}.dim == 1
            plot(levelSets{i},1,'r');
        elseif levelSets{i}.dim == 2
            plot(levelSets{i},[1,2],'r');
        else
            plot(levelSets{i},[1,2,3],'r');
        end
        resvec(end+1) = true;
    catch
        resvec(end+1) = false;
    end

    % close figure
    close;
end

figure;

% test additional plotting arguments
ls = levelSets{2}; % with '<='

plot(ls,[1,2],'PlotMethod','outer');
plot(ls,[1,2],'PlotMethod','inner');
resvec(end+1) = true;

plot(ls,[1,2],'Splits',3);
plot(ls,[1,2],'Splits',5);
resvec(end+1) = true;

plot(ls,[1,2],'Splits',3,'PlotMethod','outer');
plot(ls,[1,2],'Splits',3,'PlotMethod','outer','FaceColor','b');
plot(ls,[1,2],'Splits',3,'FaceColor','b','PlotMethod','outer');
plot(ls,[1,2],'FaceColor','b','Splits',3,'PlotMethod','outer');
resvec(end+1) = true;

% close figure
close;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
