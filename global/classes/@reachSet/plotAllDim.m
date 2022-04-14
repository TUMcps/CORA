function hans = plotAllDim(R,varargin)
% plotAllDim - Plots all 2-dimensional projections of a reachable set
%    note: only for one branch (non-hybrid, non-splitting) until now
%
% Syntax:  
%    hans = plotAllDim(R)
%    hans = plotAllDim(R,dims)
%    hans = plotAllDim(R,dims,linespec)
%    hans = plotAllDim(R,dims,linespec,'Handle',h)
%    hans = plotAllDim(R,dims,linespec,'Order',order)
%    hans = plotAllDim(R,dims,<Name-Value-pairs>)
%
% Inputs:
%    R - reachSet object
%    dims - projected dimensions to be plotted, define by
%           'all': plots all dimensions in ascending order
%           1x2 int-vector: plots 2D split
%           cell-array of 1x2 int-vectors: plots consecutive 2D splits
%    linespec - (optional) line specifications
%    h - (optional) figure handle(s) as name-value pair
%    order - (optional) zonotope order for plotting
%    splits - (optional) number of splits for refinement of polyZonotopes
%    varargin - (optional) Name-Value pairs for plotting
%    note: the last four (h, order, splits, and other Name-Value pairs)
%           can be arbitrarily ordered, linespec must (!) come before that
%
% Outputs:
%    hans - figure handle(s)
%
% Example: 
%    R; % reachSet object
%    plotAllDim(R,'all',[],'r','Order',5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Mark Wetzlinger
% Written:      23-April-2020
% Last update:  21-June-2020 (add figure handle input)
% Last revision:05-June-2020
%               18-June-2020 (remove restriction of 8 plots at max)
%               23-June-2020 (add linespec and Name-Value pairs, polyZonotope)

%------------- BEGIN CODE --------------

% default values
dims = 'all';
hans = [];
linespec = 'b';
order = [];
splits = [];
NVpairs = {};

% input processing
if nargin >= 2
    % dims has to fulfill variable requirements
    dims = varargin{1};
    if ischar(dims) && ~strcmp(dims,'all')
        error("Dimensions incorrectly defined!");
    elseif iscell(dims)
        if ~all(cellfun(@(a) all(size(a) == [1,2]),dims))
            error("Dimensions incorrectly defined!");
        else % write cells into array
            dims = cell2mat(dims);
        end
    elseif ~ischar(dims) && isvector(dims) && ~all(size(dims) == [1,2])
        error("Dimensions incorrectly defined!");
    end
    if nargin >= 3
        % process linespec and Name-Value pairs
        allowedChars = '-:.+o*xsd^v><phrgbcmykw';
        % determine start of name-value pairs
        if all(ismember(varargin{2},allowedChars))
            linespec = varargin{2}; startNV = 3;
        else % no linespec given, only name-value pairs
            startNV = 2;
        end
        % remaining inputs are name-value pairs
        NVpairs = varargin(startNV:end);
        % read order and splits
        for i = 1:2:length(NVpairs)-1
            if ischar(NVpairs{i})
                if strcmp(NVpairs{i},'Order')
                    order = NVpairs{i+1}; NVpairs{i} = []; NVpairs{i+1} = [];
                elseif strcmp(NVpairs{i},'Splits')
                    splits = NVpairs{i+1}; NVpairs{i} = []; NVpairs{i+1} = [];
                elseif strcmp(NVpairs{i},'Handle')
                    hans = NVpairs{i+1}; NVpairs{i} = []; NVpairs{i+1} = [];
                    % check figure handles
                    for h=1:length(hans)
                        if ~(ishandle(hans(h)) && strcmp(get(hans(h),'type'),'figure'))
                            error("Not a figure handle!");
                        end
                    end
                end
            end
        end
        % delete empty cells
        NVpairs = NVpairs(~cellfun('isempty',NVpairs));
    end
end

% check which dimensions should be plotted
Zdim = dim(R.timeInterval.set{1});
if strcmp(dims,'all')
    % plot all dimensions in order
    if mod(Zdim,2) == 0
        % even dimension, no repetitions
        dims = 1:Zdim;
    else
        % odd dimension, repeat second-to-last dim for last plot
        dims = [1:Zdim-1,Zdim-1,Zdim];
    end
else
    % assert that given dims are valid
    if any(dims > Zdim) ...              % dims exceed Z dimensions
            || any(dims < 1) ...         % dims smaller than 1
            || any(mod(dims,1) ~= 0) ... % dims not natural numbers
            || mod(length(dims),2) ~= 0  % number of dims are odd
        error("Given dimensions are not valid.");
    end
end


maxPerFigure = 8;
totalPlots = length(dims)/2;
totalFigures = ceil(totalPlots/maxPerFigure);
if isempty(hans)
    clear hans; % to avoid that figure returnes a number instead of figure obj
    % generate plots unless provided
    for i=1:totalFigures
        hans(i) = figure;
    end
    nrSubplots_dims = [repelem(maxPerFigure,floor(totalPlots/maxPerFigure)),...
                            mod(totalPlots,maxPerFigure)];
else
    % check if given dims work with given h
    if length(hans) ~= totalFigures
        error("Given dims do not comply with given figure handles.");
    else
        % number of subplots in each plot in h
        nrSubplots_h = arrayfun(@(a) length(get(a(1),'Children')),hans);
        % number of subplots in each plot determined by dims
        nrSubplots_dims = [repelem(maxPerFigure,floor(totalPlots/maxPerFigure)),...
                                mod(totalPlots,maxPerFigure)];
        % handling if mod is zero
        nrSubplots_dims = nrSubplots_dims(any(nrSubplots_dims));
        if ~isequal(nrSubplots_h,nrSubplots_dims)
            error("Given dims do not comply with given figure handles.");
        end
    end
end

for iFig=1:totalFigures

    figure(hans(iFig));
    
    % default for all figures before last
    rows = 2; cols = 4;
    if iFig == totalFigures
        % only last figure could have less than maxPerFigure plots
        switch nrSubplots_dims(iFig)
            case 1
                rows = 1; cols = 1;
            case 2
                rows = 1; cols = 2;
            case 3
                rows = 1; cols = 3;
            case 4
                rows = 2; cols = 2;
            case {5,6}
                rows = 2; cols = 3;
            case {7,8}
                rows = 2; cols = 4;
        end
    end

    % plotting
    for sp=1:nrSubplots_dims(iFig)
        subplot(rows,cols,sp); hold on; box on;

        xDim = dims((iFig-1)*maxPerFigure + 2*(sp-1)+1);
        yDim = dims((iFig-1)*maxPerFigure + 2*sp);

        % loop over all reachable sets
        for i=1:length(R.timeInterval.set)

            % project zonotope
            Zproj = project(R.timeInterval.set{i},[xDim,yDim]);
            
            % order reduction
            if ~isempty(order)
                Zproj = reduce(Zproj,'girard',order);
            end

            % plot set
            if isa(Zproj,'polyZonotope')
                % polyZonotope
                if isempty(splits)
                    % plot as zonotope if no splits given (much faster)
                    plot(zonotope(Zproj),[1,2],linespec,NVpairs{:});
                else
                    plot(Zproj,[1,2],linespec,'Splits',splits,NVpairs{:});
                end
            else
                % zonotope (and others?)
                plot(Zproj,[1,2],linespec,NVpairs{:});
            end
        end

        % write axis labels
        xlabel(['x_{',num2str(xDim),'}']);
        ylabel(['x_{',num2str(yDim),'}']);

    end

end


end

%------------- END OF CODE --------------