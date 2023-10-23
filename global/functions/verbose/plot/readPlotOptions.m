function NVpairs = readPlotOptions(plotOptions,varargin)
% readPlotOptions - reads all plot options (LineSpecification and
%    name-value pairs) and returns only non-redundant name-value pairs
%    Colors can be numeric, hex, matlab colors (e.g., 'r') or 'next' for
%    next color in colororder.
%
% Syntax:
%    NVpairs = readPlotOptions(plotOptions)
%
% Inputs:
%    plotOptions - LineSpecification options + Name-Value pairs
%    purpose - information about plot, admissible values:
%                   'simResult' - simulation results
%                   'reachSet' - reachability results
%                   'initialSet' - reachability results
%                   'polygon' - internal polygon class
%                   'spec:safeSet' - specification of type 'safeSet'
%                   'spec:unsafeSet' - specification of type 'unsafeSet'
%                   'spec:invariant' - specification of type 'invariant'
%                   'fill' - patches
%                   'contour' - contour
%                   'mesh' - mesh plots
%                   'surf' - surf plots
%                   'none' (default) - usual proceedings
%
% Outputs:
%    NVpairs - Name-Value pairs for Matlab plots
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/plot (all classes), colororder

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       14-July-2020 
% Last update:   29-October-2021 (only name-value pairs as output args)
%                01-June-2022 (allow empty plotOptions as input argument)
%                28-February-2023 (TL, use axis default colors)
%                24-March-2023 (TL, defaultPlotColor)
%                04-October-2023 (TL, aux_correctColorToNumeric)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input argument validation
purpose = setDefaultValues({'none'},varargin);

% check input arguments
inputArgsCheck({{purpose,'str',{'simResult','reachSet', 'initialSet', ...
        'polygon', 'spec:safeSet','spec:unsafeSet', 'spec:invariant', ...
        'mesh','surf','contour','fill','none'}}});


% allowed colors
allowed.Color = {'b','g','r','c','m','y','k','w'};

% allowed markers
allowed.Marker = {'o','+','*','.','x','_','|','s','d','^','v','<','>','p','h'};

% allowed line styles
allowed.LineStyle = {'--','-.',':','-'};

% categories
cats = {'Color','Marker','LineStyle'};

% init linespec
linespec.Color = '';
linespec.Marker = '';
linespec.LineStyle = '';

if ~isempty(plotOptions)
    % check whether first entry in plotOptions is linespec
    firstentry = lower(plotOptions{1});
    idxfound = false;
    
    % init that linespec given
    idxNVstart = 1;
    
    % loop over all categories
    for i=1:length(cats)
        instancesfound = 0;
        cat = allowed.(cats{i});
        % loop over each entry
        for s=1:length(cat)
            if contains(firstentry,cat{s})
                % instances found in this category
                instancesfound = instancesfound + 1;
                % remove found entry
                firstentry = erase(firstentry,cat{s});
                % break if two instances found in one category
                if instancesfound > 1
                    idxfound = true; break;
                end
                % overwrite empty linespec entry
                linespec.(cats{i}) = cat{s};
            end
        end
        if isempty(firstentry)
            idxfound = true; idxNVstart = 2;
        end
        if idxfound; break; end
    end
    
    % re-init linespec
    if ~isempty(firstentry)
        linespec.Color = '';
        linespec.Marker = '';
        linespec.LineStyle = '';
    end
    
    % if idxNVstart = 2, check which ones are given and overwrite by
    % corresponding NVpairs if they are provided
    plotOptions = plotOptions(idxNVstart:end);

end

% init NVpairs
NVpairs = {};
% loop over linespec
for i=1:length(cats)
    cat = cats{i};
    linespecGiven = ~isempty(linespec.(cat));
    if linespecGiven
        NVpairs = [NVpairs, cat, linespec.(cat)];
    end
    % if plotOptions given, this has precedence over linespec (if both)
    for j=1:2:length(plotOptions)
        if ischar(plotOptions{j}) && strcmp(plotOptions{j},cat)
            % assign value of name-value pair
            if linespecGiven
                NVpairs{end} = plotOptions{j+1};
            else
                NVpairs = [NVpairs, plotOptions(j:j+1)];
            end
            plotOptions(j:j+1) = [];
            break;
        end
    end
end

% distribute 'Color' to 'EdgeColor' and 'FaceColor', overwrite if necessary
[plotOptions,edgecolor] = readNameValuePair(plotOptions,'EdgeColor');
[plotOptions,facecolor] = readNameValuePair(plotOptions,'FaceColor');
[~,facealpha] = readNameValuePair(plotOptions,'FaceAlpha');
[plotOptions,filled] = readNameValuePair(plotOptions,'Filled');
[NVpairs,color] = readNameValuePair(NVpairs,'Color');

% print warning if filled given
if ~isempty(filled)
    warning("Name-value pair 'Filled'-true|false is deprecated.");
end
% note: filled is overruled if 'FaceColor' is provided
if strcmp(facecolor, 'default')
    facecolor = defaultPlotColor();
end

% set FaceColor if FaceAlpha is provided
if ~isempty(facealpha)
    % facecolor order: 'FaceColor' > 'Color' > 'EdgeColor' > default color
    if ~isempty(facecolor)
        % facecolor already set
    elseif ~isempty(color)
        facecolor = color;
    elseif ~isempty(edgecolor)
        facecolor = edgecolor;
    else
        facecolor = defaultPlotColor();
    end
end


% different handling depending on object that will be plotted
switch purpose

    case 'reachSet'
        % facecolor order: 'FaceColor' > 'Color' > default color
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'FaceColor', color];
        elseif (~isempty(filled) && ~filled) || ~isempty(edgecolor)
            NVpairs = [NVpairs, 'FaceColor', 'none'];
        else
            NVpairs = [NVpairs, 'FaceColor', defaultPlotColor()];
        end
        % edgecolor order: 'EdgeColor' > 'Color' > same as facecolor
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        elseif ~isempty(facecolor)
            NVpairs = [NVpairs, 'EdgeColor', facecolor];
        else
            NVpairs = [NVpairs, 'EdgeColor', defaultPlotColor()];
        end

    case 'initialSet'
        % edgecolor order: 'EdgeColor' > 'k'
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        else
            NVpairs = [NVpairs, 'EdgeColor', [0 0 0]];
        end
        
        % facecolor order: 'FaceColor' > 'Color' > default color
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'FaceColor', color];
        else
            NVpairs = [NVpairs, 'FaceColor', defaultPlotColor()];
        end
    
    case 'simResult'
        % just color
        if ~isempty(color)
            NVpairs = [NVpairs, 'Color', color];
        else
            NVpairs = [NVpairs, 'Color', defaultPlotColor()];
        end
   
    case 'polygon'        
        % facecolor order: 'FaceColor' > 'none'
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
        else
            NVpairs = [NVpairs, 'FaceColor', 'none'];
        end
        % edgecolor order: 'EdgeColor' > 'Color' > same as facecolor
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        elseif ~isempty(facecolor)
            NVpairs = [NVpairs, 'EdgeColor', facecolor];
        else
            NVpairs = [NVpairs, 'EdgeColor', defaultPlotColor()];
        end
        
        % Marker not supported by polygon
        [NVpairs,marker] = readNameValuePair(NVpairs,'Marker');
        if ~isempty(marker)
            warning("Name-value pair 'Marker'-<options> is ignored.");
        end

    case {'spec:safeSet','spec:invariant'}
        
        % facecolor order: 'FaceColor' > 'Color' > CORAcolor('CORA:safe')
        if ~isempty(facecolor)
            safeColor = facecolor;
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            safeColor = color;
        else
            % 'CORA:safe' same color as 'CORA:invariant'
            safeColor = CORAcolor('CORA:safe');
        end
        NVpairs = [NVpairs, 'FaceColor', safeColor];

        % edgecolor order: 'EdgeColor' > 'Color' > safeColor
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        else
            NVpairs = [NVpairs, 'EdgeColor', safeColor];
        end

    case 'spec:unsafeSet'
        
        % facecolor order: 'FaceColor' > 'Color' > CORAcolor('CORA:unsafe')
        if ~isempty(facecolor)
            unsafeColor = facecolor;
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            unsafeColor = color;
        else
            unsafeColor = CORAcolor('CORA:unsafe');
        end
        NVpairs = [NVpairs, 'FaceColor', unsafeColor];

        % edgecolor order: 'EdgeColor' > 'Color' > unsafeColor
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        else
            NVpairs = [NVpairs, 'EdgeColor', unsafeColor];
        end

    case 'fill'
    
        % facecolor order: 'FaceColor' > 'Color' > default'
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'FaceColor', color];
        else
            NVpairs = [NVpairs, 'FaceColor', defaultPlotColor()];
        end
        % edgecolor order: 'EdgeColor' > 'none'
        if isempty(edgecolor)
            NVpairs = ['EdgeColor', 'none', NVpairs];
        else
            NVpairs = ['EdgeColor', edgecolor, NVpairs];
        end
    
    case {'contour','mesh','surf'}
    
        % facecolor order: not used
        % edgecolor order: 'EdgeColor' > 'Color' > default
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        else
            NVpairs = [NVpairs, 'EdgeColor', defaultPlotColor()];
        end

    case 'none'

        % return either 'FaceColor','EdgeColor' OR 'Color'

        % facecolor order: 'FaceColor' > 'none'
        if ~isempty(facecolor) && ~strcmp(facecolor,'none')
            NVpairs = [NVpairs, 'FaceColor', facecolor];
            % edgecolor order: 'EdgeColor' > 'Color' > same as facecolor
            if ~isempty(edgecolor)
                NVpairs = [NVpairs, 'EdgeColor', edgecolor];
            elseif ~isempty(color)
                NVpairs = [NVpairs, 'EdgeColor', color];
            else % facecolor is non-empty
                NVpairs = [NVpairs, 'EdgeColor', facecolor];
            end
            % warning if filled is false
            aux_filledWarning(filled,facecolor);
        else
            % only color: 'EdgeColor' > 'Color' > default color
            usedColor = defaultPlotColor();
            if ~isempty(edgecolor)
                usedColor = edgecolor;
            elseif ~isempty(color)
                usedColor = color;
            end
            % 'Filled', true... but no facecolor -> use same as for edge
            if ~isempty(filled) && filled
                NVpairs = [NVpairs, 'FaceColor', usedColor];
                NVpairs = [NVpairs, 'EdgeColor', usedColor];
            else
                NVpairs = [NVpairs, 'Color', usedColor];
            end
        end
end

% convert char hex color to dec due to 'fill'
NVpairs = aux_correctColorToNumeric(NVpairs, 'FaceColor');
NVpairs = aux_correctColorToNumeric(NVpairs, 'EdgeColor');
NVpairs = aux_correctColorToNumeric(NVpairs, 'Color');

% add remaining plotOptions
NVpairs = [NVpairs, plotOptions];

end


% Auxiliary functions -----------------------------------------------------

function aux_filledWarning(filled,facecolor)
% print warning that name-value pair 'Filled',false is overwritten if the
% name-value pair 'FaceColor'-<color> is given (unless <color>='none')

if isempty(filled)
    return;
end

if ~filled && ~isempty(facecolor)
    if (ischar(facecolor) && ~strcmp(facecolor,'none')) || isnumeric(facecolor)
        warning("Name-value pair 'Filled'-false is overwritten by 'FaceColor'.");
    end
end

end

function NVpairs = aux_correctColorToNumeric(NVpairs, label)
    [NVpairs,color] = readNameValuePair(NVpairs,label);
    
    if~isempty(color)
        % try to convert char to numeric
        if ischar(color) || isstring(color)
            % keep old color if unsuccessful
            % error might be thrown during Matlab plot/fill
            color_new = char(color);

            % check if starts with '#'
            if color_new(1) == '#'
                 % remove '#'
                color_new = color_new(2:end);

                if length(color_new) == 6 || length(color_new) == 3
                    % check if 'abc' or 'abcdef'
                    color_new = reshape(color_new', [], 3)';
                    if size(color_new, 2) == 1
                        % extend 'abc' to 'aabbcc'
                        color_new = [color_new, color_new];
                    end
    
                    % convert hex to dec
                    color = hex2dec(color_new)'/255;
                end
            else
                % convert single character to numeric
                % this conversion is not necessary, but can cause weird
                % behavior with the current color index
                switch color
                    case 'red'
                        color = [1 0 0];
                    case 'r'
                        color = [1 0 0];
                    case 'green'
                        color = [0 1 0];
                    case 'g'
                        color = [0 1 0];
                    case 'blue'
                        color = [0 0 1];
                    case 'b'
                        color = [0 0 1];
                    case 'cyan'
                        color = [0 1 1];
                    case 'c'
                        color = [0 1 1];
                    case 'magenta'
                        color = [1 0 1];
                    case 'm'
                        color = [1 0 1];
                    case 'yellow'
                        color = [1 1 0];
                    case 'y'
                        color = [1 1 0];
                    case 'black'
                        color = [0 0 0];
                    case 'k'
                        color = [0 0 0];
                    case 'white'
                        color = [1 1 1];
                    case 'w'
                        color = [1 1 1];
                    case 'next'
                        color = defaultPlotColor();
                    otherwise
                        % keep as is
                end
            end
        end

        % add back to NVpairs
        NVpairs = [NVpairs, label, color];
    end
end

% ------------------------------ END OF CODE ------------------------------
