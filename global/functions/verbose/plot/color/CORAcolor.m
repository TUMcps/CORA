function color = CORAcolor(identifier, varargin)
% CORAcolor - returns the CORA default colors by identifier
%
% Syntax:
%    color = CORAcolor(i)
%    color = CORAcolor(identifier)
%    color = CORAcolor(identifier,varargin)
%
% Inputs:
%    i - numeric, maps to CORA:color<i>
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
%       - 'CORA:color<i>': matlab default color order colors
%       - 'CORA:blue','CORA:red',...: matlab default colors 
%       - 'MATLAB:color<i>': matlab default color order colors
%       - 'MATLAB:blue','MATLAB:red',...: matlab default colors 
%       - '...:light' for light variant 
%       - '...:dark' for dark variant
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

% Authors:       Tobias Ladner, Lukas Koller
% Written:       02-March-2023
% Last update:   24-March-2023 (TL, added 'CORA:next')
%                25-June-2024 (TL, added 'CORA:color')
%                25-August-2025 (LK, added '...:light|dark')
%                18-September-2025 (TL, enabled CORAcolor(i) as shorthand)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,Inf);

% quick exit for numeric input
if isnumeric(identifier)
    color = CORAcolor(sprintf('CORA:color%i',identifier));
    return
end

% full check 
admissableColors = {
    ... % CORA special colors
    'CORA:reachSet', 'CORA:initialSet', 'CORA:finalSet', 'CORA:simulations', ...
    'CORA:unsafe','CORA:unsafeLight','CORA:safe', 'CORA:invariant', ...
    'CORA:highlight1','CORA:highlight2', 'CORA:next', ...
    ... % default CORA colors
    'CORA:color1','CORA:blue','CORA:color2','CORA:red', ...
    'CORA:color3','CORA:yellow','CORA:color4','CORA:purple', ...
    'CORA:color5','CORA:green','CORA:color6','CORA:light-blue', ...
    'CORA:color7','CORA:dark-red', ...
    ... % default MATLAB colors
    'MATLAB:color1','MATLAB:blue','MATLAB:color2','MATLAB:red', ...
    'MATLAB:color3','MATLAB:yellow','MATLAB:color4','MATLAB:purple', ...
    'MATLAB:color5','MATLAB:green','MATLAB:color6','MATLAB:light-blue', ...
    'MATLAB:color7','MATLAB:dark-red', ...
    ... % postfixes
    '...:light','...:dark'};
inputArgsCheck({{identifier,'str',admissableColors}})

% Split the identifiers.
[pre,name,post] = aux_splitIdentifiers(identifier);

% Initialize the default color.
color = [0 0 0];

switch [pre ':' name]
    case 'CORA:reachSet'
        % varargin - {numColors, cidx}
        %    - numColors: number of reachSet colors
        %    - cidx: index of reachSet color
        narginchk(1,3);
        
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

    % cora default colors 
    % (different to matlab default colors since R2025a)
    case {'CORA:color1','CORA:blue'}                                 % blue
        color = [0 0.4470 0.7410];     
    case {'CORA:color2','CORA:red'}                                   % red
        color = [ 0.8500 0.3250 0.0980]; 
    case {'CORA:color3','CORA:yellow'}                             % yellow
        color = [0.9290 0.6940 0.1250];  
    case {'CORA:color4','CORA:purple'}                             % purple
        color = [0.4940 0.1840 0.5560];  
    case {'CORA:color5','CORA:green'}                               % green
        color = [0.4660 0.6740 0.1880];  
    case {'CORA:color6','CORA:light-blue'}                     % light blue
        color = [0.3010 0.7450 0.9330];  
    case {'CORA:color7','CORA:dark-red'}                         % dark red
        color = [0.6350 0.0780 0.1840];  

    % Matlab default colors 
    case {'MATLAB:color1','MATLAB:blue'}                             % blue
        if isMATLABReleaseOlderThan('R2025a')
            color = [0 0.4470 0.7410];
        else
            color = [0.0660 0.4430 0.7450];
        end
    case {'MATLAB:color2','MATLAB:red'}                               % red
        if isMATLABReleaseOlderThan('R2025a')
        color = [ 0.8500 0.3250 0.0980];
        else
            color = [0.8660 0.3290 0];
        end
    case {'MATLAB:color3','MATLAB:yellow'}                         % yellow
        color = [0.9290 0.6940 0.1250];

    case {'MATLAB:color4','MATLAB:purple'}                         % purple
        if isMATLABReleaseOlderThan('R2025a')
            color = [0.4940 0.1840 0.5560];
        else
            color = [0.5210 0.0860 0.8190];
        end
    case {'MATLAB:color5','MATLAB:green'}                           % green
        if isMATLABReleaseOlderThan('R2025a')
        color = [0.4660 0.6740 0.1880];
        else
            color = [0.2310 0.6660 0.1960];
        end
    case {'MATLAB:color6','MATLAB:light-blue'}                 % light blue
        if isMATLABReleaseOlderThan('R2025a')
            color = [0.3010 0.7450 0.9330];
        else
            color = [0.1840 0.7450 0.9370];
        end
    case {'MATLAB:color7','MATLAB:dark-red'}                     % dark red
        if isMATLABReleaseOlderThan('R2025a')
            color = [0.6350 0.0780 0.1840]; 
        else
            color = [0.8190 0.0150 0.5450];
        end

    % unknown color
    otherwise
        throw(CORAerror("CORA:wrongValue","first",admissableColors))
    
end

% apply color variant
color = colorvariant(post,color);
    
end


% Auxiliary functions -----------------------------------------------------

function [pre,name,post] = aux_splitIdentifiers(ids)
    % Use a regex to split the identifiers.
    tokens = regexp(ids,':','split');
    % Extract prefix.
    pre = tokens{1};
    % Extract the color name.
    name = tokens{2};
    % Extract postfix.
    if length(tokens) > 2
        post = tokens{3};
    else
        post = 'none';
    end
end

% ------------------------------ END OF CODE ------------------------------
