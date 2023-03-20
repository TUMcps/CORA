function han = plot(simRes,varargin)
% plot - plots a projection of the simulated trajectories
%
% Syntax:  
%    han = plot(simRes)
%    han = plot(simRes,dims)
%    han = plot(simRes,dims,type)
%
% Inputs:
%    obj - simResult object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%          'Height', <height> height of z-coordinate
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
% See also: simResult, simulateRandom

% Author:       Niklas Kochdumper, Matthias Althoff
% Written:      06-June-2020
% Last update:  28-July-2020 (MA, 3D plots added)
% Last revision:---

%------------- BEGIN CODE --------------

% set default input arguments
dims = setDefaultValues({[1,2]},varargin);

% check input arguments
inputArgsCheck({{simRes,'att',{'simResult'},{''}};
                {dims,'att',{'numeric'},{'positive','vector','integer'}}});

% parse plotting options
NVpairs = readPlotOptions(varargin(2:end),'simResult');
[NVpairs,height] = readNameValuePair(NVpairs,'Height','isscalar');
[NVpairs,whichtraj] = readNameValuePair(NVpairs,'Traj','ischar','x');

% check redundancy

% check dimension
if length(dims) < 2
    throw(CORAerror('CORA:plotProperties',2));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
end

% check which trajectory has to be plotted
whichtraj = checkTraj(simRes,whichtraj);

% save color index
oldColorIndex = gca().ColorOrderIndex;

% loop over all simulated trajectories
hold on
for i = 1:length(simRes.(whichtraj))
    if isempty(height) % no 3D plot
        if length(dims) == 2
            han_i = plot(simRes.(whichtraj){i}(:,dims(1)),...
                simRes.(whichtraj){i}(:,dims(2)),NVpairs{:});
        else
            han_i = plot3(simRes.(whichtraj){i}(:,dims(1)),...
                simRes.(whichtraj){i}(:,dims(2)),...
                simRes.(whichtraj){i}(:,dims(3)),NVpairs{:});
        end
    else
        % z values normalized to [0,1] in other plots
        zCoordinates = height*ones(length(simRes.(whichtraj){i}(:,dims(1))),1);
        han_i = plot3(simRes.(whichtraj){i}(:,dims(1)),...
            simRes.(whichtraj){i}(:,dims(2)), zCoordinates,NVpairs{:}); 
    end

    if i == 1
        han = han_i;
        % don't display subsequent plots in legend
        NVpairs = [NVpairs, {'HandleVisibility','off'}];
    end
end

% correct color index
updateColorIndex(oldColorIndex);

% show 3 dimensions ('hold on' causes projection to first two dimensions
%   to be shown if figure was not 3D before)
if length(dims) == 3
    view(3);
end

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