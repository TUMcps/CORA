function hans = plotAllOverTime(R,varargin)
% plotAllOverTime - plots desired projections of the reachable set over time
%
% Syntax:  
%    hans = plotOverTime(R)
%    hans = plotOverTime(R,dims)
%    hans = plotOverTime(R,dims,'k')
%    hans = plotOverTime(R,dims,'k','Handle',h)
%
% Inputs:
%    R - reachSet object
%    dims - dimensions that should be projected
%    h - (optional) figure handle(s)
%    varargin - (optional) plotting preferences
%
% Outputs:
%    hans - figure handle(s)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, reachSet

% Author:       Mark Wetzlinger, Niklas Kochdumper
% Written:      18-June-2020
% Last update:  21-June-2020 (add figure handle input)
% Last revision:---

%------------- BEGIN CODE --------------

% default values
dims = 'all';
hans = [];
linespec = 'b';
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
                if strcmp(NVpairs{i},'Handle')
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
try
    dim_x = dim(R.timeInterval.set{1});
catch
    % in case no timeInterval set given
    dim_x = dim(R.timePoint.set{1});
end
if strcmp(dims,'all')
    % plot all dimensions in order up to 8 (default)
    dims = 1:dim_x;
else
    % assert that given dims are valid
    if any(dims > dim_x) ...             % dims exceed state dimensions
            || any(dims < 1) ...         % dims smaller than 1
            || any(mod(dims,1) ~= 0) ... % dims not natural numbers
            || mod(length(dims),2) ~= 0  % number of dims are odd
        error("Given dimensions are not valid.");
    end
end

maxPerFigure = 8;
totalPlots = length(dims);
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

        % take correct state dimension
        xDim = dims((iFig-1)*maxPerFigure + sp);

        % loop over branches in R
        for i = 1:size(R,1)
            
            if ~isempty(R(i,1).timeInterval)
                Rset = R(i,1).timeInterval.set;
                Rtime = R(i,1).timeInterval.time;
            else
                Rset = R(i,1).timePoint.set;
                Rtime = R(i,1).timePoint.time;
            end
            
            % loop over all R in branch
            for j = 1:length(Rset)

                % get intervals
                intX = interval(project(Rset{j},xDim));
                intT = Rtime{j};

                intTX = cartProd(intT,intX);

                % plot interval
                plot(intTX,[1,2],linespec,NVpairs{:},'Filled',true);
            end
        end

        % write axis labels
        xlabel('t');
        ylabel(['x_{',num2str(xDim),'}']);

    end

end

end

%------------- END OF CODE --------------