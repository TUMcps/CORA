function NVpairs = readPlotOptions(plotOptions,varargin)
% readPlotOptions - reads all plot options (LineSpecification and
%    name-value pairs) and returns only non-redundant name-value pairs
%
% Syntax:  
%    NVpairs = readPlotOptions(plotOptions)
%
% Inputs:
%    plotOptions - LineSpecification options + Name-Value pairs
%    purpose - information about plot, admissible values:
%                   'simResult' - simulation results
%                   'reachSet' - reachability results
%                   'polygon' - internal polygon class
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
% See also: contSet/plot (all classes)

% Author:        Mark Wetzlinger
% Written:       14-July-2020 
% Last update:   29-October-2021 (only name-value pairs as output args)
%                01-June-2022 (allow empty plotOptions as input argument)
% Last revision: ---

%------------- BEGIN CODE --------------

% input argument validation
purpose = setDefaultValues({'none'},varargin);

% check input arguments
inputArgsCheck({{purpose,'str',{'simResult','reachSet','polygon',...
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
[plotOptions,filled] = readNameValuePair(plotOptions,'Filled');
[NVpairs,color] = readNameValuePair(NVpairs,'Color');

% print warning if filled given
if ~isempty(filled)
    warning("Name-value pair 'Filled'-true|false is deprecated.");
end
% note: filled is overruled if 'FaceColor' is provided

% different handling depending on object that will be plotted
switch purpose

    case 'reachSet'
        % facecolor order: 'FaceColor' > 'Color' > default color
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
            % warning if filled is false
            filledWarning(filled,facecolor);
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'FaceColor', color];
        elseif (~isempty(filled) && ~filled) || ~isempty(edgecolor)
            NVpairs = [NVpairs, 'FaceColor', 'none'];
        else
            NVpairs = [NVpairs, 'FaceColor', colorblind('b')];
        end
        % edgecolor order: 'EdgeColor' > 'Color' > same as facecolor
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        elseif ~isempty(facecolor)
            NVpairs = [NVpairs, 'EdgeColor', facecolor];
        else
            NVpairs = [NVpairs, 'EdgeColor', colorblind('b')];
        end
    
    case 'simResult'
        % just color
        if ~isempty(color)
            NVpairs = [NVpairs, 'Color', color];
        else
            NVpairs = [NVpairs, 'Color', colorblind('y')];
        end
   
    case 'polygon'
        % only called internally, so 'Filled' not respected
        
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
            NVpairs = [NVpairs, 'EdgeColor', colorblind('b')];
        end
        
        % Marker not supported by polygon
        [NVpairs,marker] = readNameValuePair(NVpairs,'Marker');
        if ~isempty(marker)
            warning("Name-value pair 'Marker'-<options> is ignored.");
        end

    case 'fill'
    
        % facecolor order: 'FaceColor' > 'Color' > default'
        if ~isempty(facecolor)
            NVpairs = [NVpairs, 'FaceColor', facecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'FaceColor', color];
        else
            NVpairs = [NVpairs, 'FaceColor', colorblind('b')];
        end
        % edgecolor order: always 'none'
        NVpairs = [NVpairs, 'EdgeColor', 'none'];
    
    case {'contour','mesh','surf'}
    
        % facecolor order: not used
        % edgecolor order: 'EdgeColor' > 'Color' > default
        if ~isempty(edgecolor)
            NVpairs = [NVpairs, 'EdgeColor', edgecolor];
        elseif ~isempty(color)
            NVpairs = [NVpairs, 'EdgeColor', color];
        else
            NVpairs = [NVpairs, 'EdgeColor', colorblind('b')];
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
            filledWarning(filled,facecolor);
        else
            % only color: 'EdgeColor' > 'Color' > default color
            usedColor = colorblind('b'); % define like this because re-used
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

% add remaining plotOptions
NVpairs = [NVpairs, plotOptions];

end


% Auxiliary function
function filledWarning(filled,facecolor)
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

%------------- END OF CODE --------------
