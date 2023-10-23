function han = plotOverTime(specs,varargin)
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

% Authors:       Tobias Ladner
% Written:       03-March-2023
% Last update:   ---
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[specs,dims,NVpairs] = aux_parseInput(specs,varargin{:});

% 2. preprocess
[spectime,holdStatus] = aux_preprocess(specs);

% 3. plot specifications
han = aux_plotSpecs(specs,dims,NVpairs,spectime);

% 4 postprocess
aux_postprocess(holdStatus)

% 5. clear han
if nargout == 0
    clear han
end

end


% Auxiliary functions -----------------------------------------------------

function [specs,dims,NVpairs] = aux_parseInput(specs,varargin)
    % parse input

    % default values for the optional input arguments
    dims = setDefaultValues({1},varargin);
    
    % check input arguments
    inputArgsCheck({{specs,'att','specification','nonempty'};
                    {dims,'att','numeric',{'nonempty','scalar','integer','positive'}}});

    % plot options are read depending on the specification type
    NVpairs = varargin(2:end);
end

function [spectime,holdStatus] = aux_preprocess(specs)
    % preprocess

    % check if any of the specifications have a defined time frame
    spectime = cell(length(specs),1);
    for i=1:length(specs)
        spectime{i} = specs(i).time;
    end
    noTimeGiven = cellfun(@(x) representsa_(x,'emptySet',eps), ...
        spectime,'UniformOutput',true);
    if all(noTimeGiven)
        throw(CORAerror('CORA:specialError','No specification has a defined time frame.'));
    elseif any(noTimeGiven)
        % set empty .time fields to min/max of others
        tmin = Inf; tmax = -Inf;
        for i=1:length(specs)
            if ~noTimeGiven(i)
                tmin = min([tmin,spectime{i}.inf]);
                tmax = max([tmax,spectime{i}.sup]);
            end
        end
        spectime{noTimeGiven} = interval(tmin,tmax);
    end

    % check hold status
    holdStatus = ishold;
    if ~holdStatus
        % clear any plot if hold status is 'off'
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        ax = gca();
        set(ax,'ColorOrderIndex',1);
    end
    hold on;
end

function han = aux_plotSpecs(specs,dims,NVpairs,spectime)
    % plot specifications

    % reshape specs to column vector
    specs = reshape(specs,[],1);
    types = arrayfun(@(spec) string(spec.type), specs);

    % gather unique types
    uniqueTypes = unique(types,'rows')'; 
    
    for type = uniqueTypes
        % read all specs of current type
        idx = strcmp(types,type);
        specs_type = specs(idx);
        spectime_type = spectime(idx); 
        
        % read plotting options depending on type
        switch type
            case {'safeSet','unsafeSet','invariant'}
                % read specific specification plot options for these types
                NVpairs_type = readPlotOptions(NVpairs, sprintf("spec:%s",type));
    
            case 'custom'
                % read default plot options
                NVpairs_type = readPlotOptions(NVpairs);
    
            otherwise
                throw(CORAerror('CORA:notSupported',...
                sprintf("Projection a specifications of type '%s' is not yet supported.",type)));
        end

        % init
        sets = cell(1, length(specs_type));
        
        % iterate over specs of current type
        for i=1:length(sets)
            sets{i} = cartProd( ...
                spectime_type{i}, ...                           % x
                interval(project(specs_type(i).set,dims)) ...   % y
            );
        end
        
        % plot all sets
        han = plotMultipleSetsAsOne(sets,[1,2],NVpairs_type);
    end    
end

function aux_postprocess(holdStatus)
    % reset hold status
    if ~holdStatus
        hold("off");
    end
end

% ------------------------------ END OF CODE ------------------------------
