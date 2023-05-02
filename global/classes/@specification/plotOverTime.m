function han = plotOverTime(spec,varargin)
% plotOverTime - plots a projection of the specification
%
% Syntax:  
%    han = plotOverTime(spec)
%    han = plotOverTime(spec,dims)
%    han = plotOverTime(spec,dims,plotOptions)
%
% Inputs:
%    spec - specification object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Author:       Tobias Ladner
% Written:      03-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{spec,'att','specification','nonempty'};
                {dims,'att','numeric',{'nonempty','scalar','integer','positive'}}});

% check hold status
holdStatus = ishold;
if ~holdStatus
    plot(NaN,NaN,'HandleVisibility','off');
    % reset color index (before readPlotOptions!)
    set(gca(),'ColorOrderIndex',1);
end

% check if any of the specifications have a defined time frame
spectime = cell(length(spec),1);
for i=1:length(spec)
    spectime{i} = spec(i).time;
end
noTimeGiven = cellfun(@(x) isempty(x),spectime,'UniformOutput',true);
if all(noTimeGiven)
    throw(CORAerror('CORA:specialError','No specification has a defined time frame.'));
elseif any(noTimeGiven)
    % set empty .time fields to min/max of others
    tmin = Inf; tmax = -Inf;
    for i=1:length(spec)
        if ~noTimeGiven(i)
            tmin = min([tmin,spectime{i}.inf]);
            tmax = max([tmax,spectime{i}.sup]);
        end
    end
    spectime{noTimeGiven} = interval(tmin,tmax);
end

for i=1:length(spec)
    hold on;
    spec_i = spec(i);

    % read plotting options depending on type
    switch spec_i.type
        case {'safeSet','unsafeSet'}
            % check name-value pairs
            NVpairs = readPlotOptions(varargin(2:end), ...
                sprintf("spec:%s",spec_i.type));

        case 'custom'
            NVpairs = readPlotOptions(varargin(2:end));

        otherwise
            throw(CORAerror('CORA:notSupported',...
            sprintf("Projection a specifications of type '%s' is not yet supported.",spec_i.type)));
    end
    if i>1
        % don't show in legend
        NVpairs = [NVpairs,'HandleVisibility','off'];
    end

    % combine set with time
    set_i = cartProd( ...
        spectime{i}, ...                        % x
        interval(project(spec_i.set,dims)) ...  % y
    );

    % plot set
    han_i = plot(set_i, [1,2], NVpairs{:});
    if i==1
        han = han_i;
    end
end

% reset hold status
if ~holdStatus
    hold off
end

if nargout == 0
    clear han
end

end

%------------- END OF CODE --------------