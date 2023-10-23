function fig = fig2tikz(varargin)
% fig2tikz - re-plot objects in figure so that matlab2tikz can be executed;
%    currently, only polygon object are supported for conversion;
%    additionally, line objects are simplified
%
% Syntax:
%    fig = fig2tikz
%    fig = fig2tikz(fig)
%    fig = fig2tikz(res)
%    fig = fig2tikz(fig,res)
%
% Inputs:
%    fig - handle to figure object
%    res - precision of simplified objects
%
% Outputs:
%    fig - handle to figure object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-October-2021
% Last update:   03-November-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
fig = gcf;
precision = 1e-5;

% input preprocessing
if nargin > 0
    % indices
    maxIdx = length(varargin); idx = 1;
    
    % order: figure, res
    while true
        % check figure
        if isgraphics(varargin{idx}) && strcmp(varargin{idx}.Type,'figure')
            fig = varargin{idx}; idx = idx + 1;
        end
        if idx > maxIdx; break; end
        % check resolution
        if isnumeric(varargin{idx}) && isscalar(varargin{idx})
            precision = varargin{idx}; idx = idx + 1;
        end
        if idx > maxIdx; break; end
    end
end


% loop over all subplots
nrChildren = length(fig.Children);

for c=1:nrChildren

    % skip child if not an Axes object
    if ~strcmp(get(fig.Children(c),'type'),'axes')
        continue;
    end
    % select current axes
    set(gcf,'CurrentAxes',fig.Children(c));
    
    % read out number of all plotted elements
    numElem = length(fig.Children(c).Children);

    % loop over all elements in subplot (Axes object)
    for i=1:numElem
        
        % read out all elements
        elem = fig.Children(c).Children(i);
        
        % check which elements are not compatible with matlab2tikz
        switch elem.Type

            case {'patch','image','hggroup','light','', ...
                    'matlab.graphics.primitive.Group','scatter','bar',...
                    'stair','stem','errorbar','area','quiver','contour',...
                    'hgtransform','surface','text','rectangle','histogram'}

                % all supported by matlab2tikz (see line ~730)
                % no changes in this function

            case 'line'

                % supported by matlab2tikz, but we want to reduce the
                % number of data points to limit the resulting file size

                % reduced number of points
                V_red = reducepoly([elem.XData; elem.YData]',precision);

                % update XData and YData using simplified representation
                elem.XData = V_red(:,1);
                elem.YData = V_red(:,2);

            case 'polygon'

                % fix for display name (cannot be empty)
                displayname = elem.DisplayName;
                if isempty(displayname)
                    % dummy name to circumvent error in matlab2tikz
                    displayname = ['Object ' num2str(i)];
                end

                % read figure properties (incl. legend entry!)
                properties = {'FaceColor',elem.FaceColor,...
                    'FaceAlpha',elem.FaceAlpha,...
                    'EdgeColor',elem.EdgeColor,...
                    'LineStyle',elem.LineStyle,...
                    'LineWidth',elem.LineWidth,...
                    'DisplayName',displayname};

                % read vertices
                V = elem.Shape.Vertices;

                % NaN entries in V serve as a delimiter for separate lines
                idxNaN = find(isnan(V(:,1)),1);

                if ~isempty(idxNaN)
                    % currently only one NaN entry supported which corresponds
                    % to a ring (one closed line inside another)
                    if length(idxNaN) == 1
                        % start part and end part
                        V1 = V(1:idxNaN-1,:);
                        V2 = V(idxNaN+1:end,:);

                        % find out which line is which
                        if all([min(V1,[],1) < min(V2,[],1), max(V1,[],1) > max(V2,[],1)])
                            % order is ok
                        elseif all([min(V1,[],1) > min(V2,[],1), max(V1,[],1) < max(V2,[],1)])
                            % swap order
                            temp = V1; V1 = V2; V2 = temp;
                        else
                            % unknown shape (not a ring)
                            throw(CORAerror('CORA:specialError','Unknown shape'));
                        end

                        % reduce resolution of objects
                        V_red1 = reducepoly(V1,precision);
                        V_red2 = reducepoly(V2,precision);

                        % set 'FaceColor' of inner plot to white
                        properties2 = properties;
                        properties2{2} = 'w';

                        % plot matlab2tikz-compatible objects
                        patch('XData',V_red1(:,1),'YData',V_red1(:,2),properties{:});
                        patch('XData',V_red2(:,1),'YData',V_red2(:,2),properties2{:});

                        % delete previous graphics object
                        delete(elem);

                        % remove legend entry of inner plot: this is a
                        % workaround caused by MATLAB as one cannot plot a
                        % patch without generating an empty legend entry
                        % (e.g., by omitting 'DisplayName' or setting it to '')
                        fig.Children(c).Legend.String = ...
                            fig.Children(c).Legend.String(1:end-1);

                        % note: new children are appended to the beginning of the list,
                        % thus we require to reorder the plots after replotting
                        fig.Children(c).Children = ...
                            fig.Children(c).Children([3:i+1 1:2 i+2:numElem+1]);

                        % update number of elements
                        numElem = length(fig.Children(c).Children);

                    else
                        throw(CORAerror('CORA:notSupported',...
                            'This type of polygon is currently not supported.'));
                    end

                else
                    % reduce resolution of object
                    V_red = reducepoly(V,precision);

                    % plot matlab2tikz-compatible object
                    patch('XData',V_red(:,1),'YData',V_red(:,2),properties{:});

                    % delete previous graphics object
                    delete(elem);

                    % note: new children are appended to the beginning of the list,
                    % thus we require to reorder the plots after replotting
                    fig.Children(c).Children = ...
                        fig.Children(c).Children([2:i 1 i+1:numElem]);
                end

            otherwise

                % not supported
                throw(CORAerror('CORA:notSupported',...
                    ['Object of class ' elem.Type ' are currently not supported.']));

        end

    end

end

% scale precision of all elements
% cleanfigure('scalePrecision',precision);

% return figure
fig = gcf;

% ------------------------------ END OF CODE ------------------------------
