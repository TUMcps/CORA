function res = testLongDuration_levelSet_plot
% testLongDuration_levelSet_plot - unit test function of plot
%
% Syntax:  
%    res = testLongDuration_levelSet_plot
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

% Author:       Maximilian Perschl
% Written:      08-November-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

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
    catch
        close;
        res = false;
    end

end

close all;

%------------- END OF CODE --------------