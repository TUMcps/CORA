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
% Last update:   24-March-2023 (TL, 'CORA:next')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 1
    throw(CORAerror("CORA:notEnoughInputArgs", 1))
end
inputArgsCheck({{identifier, 'str', {'CORA:reachSet', ...
    'CORA:initialSet', 'CORA:finalSet', 'CORA:simulations', ...
    'CORA:unsafe','CORA:safe', 'CORA:invariant', ...
    'CORA:highlight1','CORA:highlight2', 'CORA:next'}}})

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
    case {'CORA:safe','CORA:invariant'}
        color = [0.4706 0.7725 0.4980]; % green
    case 'CORA:highlight1'
        color = [1.0000 0.6824 0.2980]; % orange
    case 'CORA:highlight2'
        color = [0.6235 0.7294 0.2118]; % light green
    case 'CORA:next'
        color = defaultPlotColor();
        
end
    
end

% ------------------------------ END OF CODE ------------------------------
