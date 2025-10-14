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

% Authors:       Maximilian Perschl, Niklas Kochdumper
% Written:       08-November-2021
% Last update:   09-September-2025 (NK, added random case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 0.01;

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
% 3d multiple equality
eq = [40*x*y^3*z^3;-56*x*z^6 - 48*x*y^2*z^3 - 180*x^2*y^3*z^2];
levelSets{6} = levelSet(eq,[x;y;z],'==');
% 1d plot:
eq = x^2 - 4;
levelSets{7} = levelSet(eq,x,'==');
% empty set
eq = sym(7);
levelSets{8} = levelSet(eq,x,'==');

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
    catch ME
        rethrow(ME)
    end

    % close figure
    close;
end

figure;

% test additional plotting arguments
ls = levelSets{2}; % with '<='

plot(ls,[1,2],'PlotMethod','outer');
plot(ls,[1,2],'PlotMethod','inner');

plot(ls,[1,2],'Splits',3);
plot(ls,[1,2],'Splits',5);

plot(ls,[1,2],'Splits',3,'PlotMethod','outer');
plot(ls,[1,2],'Splits',3,'PlotMethod','outer','FaceColor','b');
plot(ls,[1,2],'Splits',3,'FaceColor','b','PlotMethod','outer');
plot(ls,[1,2],'FaceColor','b','Splits',3,'PlotMethod','outer');

% close figure
close;

% 2D Inequality Constraints (random) ----------------------------------

% generate random point cloud
N = 10 + randi(20);
s = randi(10);

points = s*(-1 + 2*rand(2,N)) + s*(-1 + 2*rand(2,1));

% generate level set enclosing the point cloud
method = {'single','multiple'};
order = 4;

ls = levelSet.enclosePoints(points,method{randi(2)},order);

% plot level set
figure; hold on;
xlim([min(points(1,:))-1,max(points(1,:)+1)]); 
ylim([min(points(2,:))-1,max(points(2,:)+1)]);
h = plot(ls);

% extract filled area from the plot
pgon = [];

for i = 1:length(h)
    pgon = pgon | polygon(h(i).Vertices(:,1),h(i).Vertices(:,2));
end

% generate random points inside and outside the level set
pIn = randPoint(pgon,100);
pOut = randPoint(subtract(polygon(interval(pgon)),pgon),100);

minDist = max(rad(interval(pgon)))/20;

indIn = []; indOut = [];

for i = 1:size(pIn,2)
    if hausdorffDist(pgon,pIn(:,i)) > minDist
        indIn = [indIn,i];
    end
    if hausdorffDist(pgon,pOut(:,i)) > minDist
        indOut = [indOut,i];
    end
end

pIn = pIn(:,indIn); pOut = pOut(:,indOut);

% check for correctness
assert(all(contains(ls,pIn,'exact',tol)));

if ~isempty(pOut)
    if length(ls.eq) == 1
        assert(all(contains(not(ls),pOut,'exact',tol)));
    else
        assert(all(~contains(ls,pOut)));
    end
end

close;

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
