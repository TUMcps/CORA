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
inputArgsCheck({{spec,'att',{'specification'},{''}};
                {dims,'att',{'numeric'},{'nonempty','scalar','integer','positive'}}});

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
    set = cartProd( ...
        spec_i.time, ...                        % x
        interval(project(spec_i.set,dims)) ...  % y
    );

    % plot set
    han_i = plot(set, [1,2], NVpairs{:});
    if i==1
        han = han_i;
    end
end

if nargout == 0
    clear han
end

end

%------------- END OF CODE --------------