function han = plotOverTime(simRes,varargin)
% plotOverTime - plots the simulated trajectories over time
%
% Syntax:  
%    han = plotOverTime(simRes)
%    han = plotOverTime(simRes,dims)
%    han = plotOverTime(simRes,dims,type)
%
% Inputs:
%    simRes - simResult object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%          'Traj', <whichtraj> corresponding to
%                   x ... state trajectory (default)
%                   y ... output trajectory (default if no ti)
%                   a ... algebraic trajectory
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: simResult, simulateRandom, plot

% Author:       Niklas Kochdumper
% Written:      06-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = setDefaultValues({1},varargin);

% parse input arguments
NVpairs = readPlotOptions(varargin(2:end),'simResult');
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% check which trajectory has to be plotted
whichtraj = checkTraj(simRes,whichtraj);

% save color index
oldColorIndex = gca().ColorOrderIndex;

% loop over all simulated trajectories
hold on
for i = 1:length(simRes.(whichtraj))
    han_i = plot(simRes.t{i},simRes.(whichtraj){i}(:,dims),NVpairs{:});

    if i == 1
        han = han_i;
        % don't display subsequent plots in legend
        NVpairs = [NVpairs, {'HandleVisibility','off'}];
    end
end

% correct color index
updateColorIndex(oldColorIndex);

if nargout == 0
    clear han;
end

end


% Auxiliary function ------------------------------------------------------

function whichtraj = checkTraj(simRes,whichtraj)

% must be character vector for switch-expression to work properly
if isempty(whichtraj)
    % default value
    whichtraj = 'x';
end

switch whichtraj
    case 'x'
        % no issues (should always be there)
        
    case 'y'
        if isempty(simRes(1).y)
            throw(CORAerror('CORA:emptyProperty'));
        end

    case 'a'
        if isempty(simRes(1).a)
            throw(CORAerror('CORA:emptyProperty'));
        end

    otherwise
        % error
        throw(CORAerror('CORA:specialError','Wrong value for name-value pair'));
end

end

%------------- END OF CODE --------------