function color = CORAcolor(identifier, varargin)
% CORAcolor - returns the CORA default colors by identifier
%
% Syntax:
%    color = CORAcolor(identifier)
%
% Inputs:
%    identifier - name of CORA colormap, one of:
%       - 'CORA:reachSet'
%       - 'CORA:initialSet'
%       - 'CORA:finalSet'
%       - 'CORA:simulations'
%       - 'CORA:unsafe'
%       - 'CORA:safe'
%       - 'CORA:highlight1': orange
%       - 'CORA:highlight2': blue
%       - 'CORA:next': next color according to colororder
%    varargin - depends on identifier, see below
%
% Outputs:
%    color - rbg triple
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: useCORAcolors, readPlotOptions

% Authors:       Tobias Ladner
% Written:       02-March-2023
% Last update:   24-March-2023 (TL, added 'CORA:next')
%                25-June-2024 (TL, added 'CORA:color')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 1
    throw(CORAerror("CORA:notEnoughInputArgs", 1))
end
inputArgsCheck({{identifier, 'str', {'CORA:reachSet', ...
    'CORA:initialSet', 'CORA:finalSet', 'CORA:simulations', ...
    'CORA:unsafe','CORA:unsafeLight','CORA:safe', 'CORA:invariant', ...
    'CORA:highlight1','CORA:highlight2', 'CORA:next', ...
    'CORA:color1','CORA:color2','CORA:color3','CORA:color4', ...
    'CORA:color5','CORA:color6','CORA:color7'}}})

color = [0 0 0]; % default

switch identifier
    case 'CORA:reachSet'
        % varargin - {numColors, cidx}
        %    - numColors: number of reachSet colors
        %    - cidx: index of reachSet color
        if nargin > 3
            throw(CORAerror('CORA:tooManyInputArgs', 3))
        end
        
        [numColors, cidx] = setDefaultValues({1 1}, varargin);
        inputArgsCheck({ ...
            {identifier, 'str', 'CORA:reachSet'}; ... 
            {numColors, 'att', 'numeric', {'integer','scalar','positive'}};
            {cidx, 'att', 'numeric', {'integer','scalar','positive'}};
            })
        if cidx > numColors
            throw(CORAerror('CORA:wrongValue', 'second/third', 'Color index must not be larger than number of colors.'))
        end

        colorMain  = [0.2706 0.5882 1.0000]; % blue
        colorWorse = [0.6902 0.8235 1.0000]; % light blue

        if cidx == numColors
            color = colorMain;
        elseif cidx == 1
            color = colorWorse;
        else
            color = colorWorse + (colorMain-colorWorse)*((cidx-1)/(numColors-1));
        end

    case 'CORA:initialSet'
        color = [1 1 1];
    case 'CORA:finalSet'
        color = [1 1 1] * 0.9;
    case 'CORA:simulations'
        color = [0 0 0];
    case 'CORA:unsafe'
        color = [0.9451 0.5529 0.5686]; % red
    case 'CORA:unsafeLight'
        color = [0.9059 0.7373 0.7373]; % light red
    case {'CORA:safe','CORA:invariant'}
        color = [0.4706 0.7725 0.4980]; % green
    case 'CORA:highlight1'
        color = [1.0000 0.6824 0.2980]; % orange
    case 'CORA:highlight2'
        color = [0.6235 0.7294 0.2118]; % light green
    case 'CORA:next'
        color = defaultPlotColor();

    % matlab default colors
    case 'CORA:color1'
        color = [0 0.4470 0.7410];       % blue
    case 'CORA:color2'
        color = [ 0.8500 0.3250 0.0980]; % red
    case 'CORA:color3'
        color = [0.9290 0.6940 0.1250];  % yellow
    case 'CORA:color4'
        color = [0.4940 0.1840 0.5560];  % purple
    case 'CORA:color5'
        color = [0.4660 0.6740 0.1880];  % green
    case 'CORA:color6'
        color = [0.3010 0.7450 0.9330];  % light blue
    case 'CORA:color7'
        color = [0.6350 0.0780 0.1840];  % dark red
        
end
    
end

% ------------------------------ END OF CODE ------------------------------
