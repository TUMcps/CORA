function han = plotMultipleSetsAsOne(sets,varargin)
% plotMultipleSetsAsOne - plots a list of sets to the current figure
%    while making it appear as a single set
%    (e.g. in legend, colororder, hold status, ...)
%
% Syntax:
%    han = plotMultipleSetsAsOne(V,varargin)
%
% Inputs:
%    sets - matrix storing the polygon vertices
%    dims - desired dimensions to plot
%    NVpairs - plot settings name-value pairs
%           or cell array containing NVpairs for each set
%           (only first entry of cell array is shown in legend)
%    purpose - purpose - information about plot
%           or cell array containing purpose for each set
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/plot, readPlotOptions

% Authors:       Tobias Ladner
% Written:       12-July-2023
% Last update:   04-October-2023 (TL, catch polyshape objects)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. prepare current axis (before parsing input)
[holdStatus,oldColorIndex] = aux_prepareAxis();

% 2. parse input
[dims,NVpairs] = aux_parseInput(sets,varargin{:});

% 3. plot sets
han = aux_plotSets(sets,dims,NVpairs);

% 4. postprocess axis
aux_postprocessAxis(holdStatus,oldColorIndex)

% clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [dims,NVpairs] = aux_parseInput(sets,varargin)
    % parse input
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs', 4))
    end

    % set default values
    [dims,NVpairs,purpose] = setDefaultValues({[1,2],{},'none'},varargin);

    % check input args
    inputArgsCheck({{sets, 'att', 'cell'},
        {dims,'att','numeric',{'nonempty','vector','integer','positive'}}})

    % extend NVpairs to match dimension of sets
    if isempty(NVpairs) || ~iscell(NVpairs{1})
        NVpairs = repmat({NVpairs},max(1,length(sets)),1);
    end

    % check if dimensions match
    if ~all(length(sets) == length(NVpairs)) && ~isempty(sets)
        throw(CORAerror('CORA:wrongInputInConstructor','Either specify one NVpairs object for all sets or a cell array of NVpairs for each set.'))
    end

    % extend purpose to match dimension of sets
    if ~iscell(purpose)
        purpose = repmat({purpose},max(1,length(sets)),1);
    end

    % check if dimensions match
    if ~all(length(sets) == length(purpose)) && ~isempty(sets)
        throw(CORAerror('CORA:wrongInputInConstructor','Either specify one purpose for all sets or a cell array of purposes for each set.'))
    end

    % iterate over all sets and read plot options
    for i=1:length(sets)
        % read NVpair and purpose for set i
        NVpairs_i = NVpairs{i};
        purpose_i = purpose{i};

        % read anyway if no NVpairs are given
        NVpairs_i = readPlotOptions(NVpairs_i, purpose_i);

        % store new NVpair for set i
        NVpairs{i} = NVpairs_i;
    end
end

function [holdStatus,oldColorIndex] = aux_prepareAxis()
    % prepare current axis
    ax = gca();

    % prepare hold status
    holdStatus = ishold;
    if ~holdStatus
        % clear any plot if hold status is 'off'
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        set(ax,'ColorOrderIndex',1);
    end
    hold on;

    % prepare index in color order
    oldColorIndex = ax.ColorOrderIndex;
end

function han = aux_plotSets(sets,dims,NVpairs)
    % plot sets

    if isempty(sets)
        % plot empty set
        han = plotPolygon(zeros(length(dims),0),NVpairs{1});

    else
        % default visibility (taken from NVpairs)
        handleVis = {};

        % iterate over all sets
        for i=1:length(sets)
            % read set i
            sets_i = sets{i};
            NVpairs_i = NVpairs{i};

            % plot set i
            if isa(sets_i,'polyshape')
                sets_i = polygon(sets_i);
            end
            han_i = plot(sets_i,dims,NVpairs_i{:}, handleVis{:});
            
                
            if i == 1    
                % force not showing subsequent plots in legend            
                handleVis = {'HandleVisibility','off'};

                % return first plotted set as handle
                han = han_i;
            end
        end
    end
end

function aux_postprocessAxis(holdStatus,oldColorIndex)
    % postprocess axis

    % correct color index
    updateColorIndex(oldColorIndex);

    % reset hold status
    if ~holdStatus
        hold("off");
    end
end

% ------------------------------ END OF CODE ------------------------------
