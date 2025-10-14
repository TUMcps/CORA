function res = test_plotMultipleSetsAsOne
% test_plotMultipleSetsAsOne - unit test function for plotMultipleSetsAsOne
%
% Syntax:
%    res = test_plotMultipleSetsAsOne
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
% See also: none

% Authors:       Tobias Ladner
% Written:       13-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create sets
S1 = zonotope([1;2],[1 0 1; 0 1 1]);
S2 = interval([-5;-3],[-3;-2]);
S3 = zonotope([-4;5],[1;0]);

try

    figure;
    ax = gca(); hold off;

    % test plot with differnt amount of sets
    plotMultipleSetsAsOne({})
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1})
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2})
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2,S3})
    assert(length(ax.Children) == 1);

    % check dimension
    plotMultipleSetsAsOne({S1,S2,S3},[1,2]);
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2,S3},[2,1]);
    assert(length(ax.Children) == 1);

    % check if all handles are returned
    han = plotMultipleSetsAsOne({S1,S2,S3},[2,1]);
    assert(numel(han) == 3);

    % check plot options
    plotMultipleSetsAsOne({S1,S2,S3},[1,2],{'FaceColor',[0 0 1]});
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2,S3},[1,2],{'FaceColor',[0 0 1],'EdgeColor',[1 0 0]});
    assert(length(ax.Children) == 1);

    % check purpose
    plotMultipleSetsAsOne({S1,S2,S3},[1,2],{},'none');
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2,S3},[1,2],{'Marker','.'},'none');
    assert(length(ax.Children) == 1);

    % check different settings for each set (required for reachSet/plot)
    plotMultipleSetsAsOne({S1,S2,S3},[1,2],{{'Marker','.'},{'Marker','o'},{'Marker','.'}});
    assert(length(ax.Children) == 1);

    close

    % test with hold on
    figure; 
    ax = gca(); hold on;

    % test plot with differnt amount of sets
    plotMultipleSetsAsOne({S1})
    assert(length(ax.Children) == 1);
    plotMultipleSetsAsOne({S1,S2})
    assert(length(ax.Children) == 2);
    plotMultipleSetsAsOne({S1,S2,S3})
    assert(length(ax.Children) == 3);

    close;

catch ME
    close
    rethrow(ME)
end

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
