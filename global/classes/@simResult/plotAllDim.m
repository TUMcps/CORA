function hans = plotAllDim(obj,varargin)
% plotAllDim - plots a desired subset of state variables over time
%    note: only for one branch (non-hybrid, non-splitting) until now
%
% Syntax:  
%    hans = plotAllDim(obj)
%    hans = plotAllDim(obj,dims)
%    hans = plotAllDim(obj,dims,'r')
%    hans = plotAllDim(obj,dims,'r','Handle',h,'LineWidth',2,...)
%
% Inputs:
%    R - reachSet object
%    dims - projected dimensions to be plotted (1x2n int-vector or 'all')
%    linespec - (optional) LineSpec properties
%    h - (optional) figure handle(s) as name-value pair
%    NVpairs - (optional) name-value pairs
%
% Outputs:
%    hans - figure handle(s)
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon

% Author:       Mark Wetzlinger
% Written:      21-June-2020
% Last update:  23-June-2020
% Last revision:---

%------------- BEGIN CODE --------------

% default values
dims = 'all';
hans = [];
linespec = 'b';
NVpairs{1} = {};

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
                if strcmp(NVpairs{i},'Handle')
                    hans = NVpairs{i+1}; NVpairs{i} = []; NVpairs{i+1} = [];
                    % check figure handles
                    for iHan=1:length(hans)
                        if ~(ishandle(hans(iHan)) && strcmp(get(hans(iHan),'type'),'figure'))
                            error("Third input must be a figure handle!");
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
dim_x = size(obj.x{1},2);
if strcmp(dims,'all')
    % plot all dimensions in order
    if mod(dim_x,2) == 0
        % even dimension, no repetitions
        dims = 1:dim_x;
    else
        % odd dimension, repeat second-to-last dim for last plot
        dims = [1:dim_x-1,dim_x-1,dim_x];
    end
else
    % assert that given dims are valid
    if any(dims > dim_x) ...              % dims exceed Z dimensions
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
    clear h; % to avoid that figure returnes a number instead of figure obj
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

        for i=1:length(obj.x)
            % plot simulation trajectory
            plot(obj.x{i}(:,xDim),obj.x{i}(:,yDim),linespec,NVpairs{:});
        end

        % write axis labels
        xlabel(['x_{',num2str(xDim),'}']);
        ylabel(['x_{',num2str(yDim),'}']);

    end

end

end

%------------- END OF CODE --------------