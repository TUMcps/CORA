function han = plot(specs,varargin)
% plot - plots a projection of the specification
%
% Syntax:
%    han = plot(spec)
%    han = plot(spec,dims)
%    han = plot(spec,dims,plotOptions)
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
holdStatus = aux_preprocess();

% 3. plot specifications
han = aux_plotSpecs(specs,dims,NVpairs);

% 4 postprocess
aux_postprocess(holdStatus)

% 5. clear han
if nargout == 0
    clear han
end

end


% Auxiliary functions -----------------------------------------------------

function [specs,dims,NVpairs] = aux_parseInput(specs,varargin)
    % default values for the optional input arguments
    dims = setDefaultValues({[1,2]},varargin);
    
    % check input arguments
    inputArgsCheck({{specs,'att','specification','nonempty'};
                    {dims,'att','numeric',{'nonempty','vector','integer','positive'}}});
    
    % check dimension
    if length(dims) < 2
        throw(CORAerror('CORA:plotProperties',2));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end

    % plot options are read depending on the specification type
    NVpairs = varargin(2:end);
end

function holdStatus = aux_preprocess()
    % preprocess

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

function han = aux_plotSpecs(specs,dims,NVpairs)
    % plot specifications, grouped per type
    
    % reshape specs to column vector
    specs = reshape(specs,[],1);
    types = arrayfun(@(spec) string(spec.type), specs);

    % gather unique types
    uniqueTypes = unique(types,'rows')'; 

    for type = uniqueTypes
        % read all specs of current type
        specs_type = specs(strcmp(types,type));
        
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
        
        % plot all sets
        han = plotMultipleSetsAsOne({specs_type.set},dims,NVpairs_type);
    end
end

function aux_postprocess(holdStatus)
    % reset hold status
    if ~holdStatus
        hold("off");
    end
end

% ------------------------------ END OF CODE ------------------------------
