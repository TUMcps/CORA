function han = plotOverTime(R,varargin)
% plotOverTime - plots the reachable set over time
%
% Syntax:
%    han = plotOverTime(R)
%    han = plotOverTime(R,dims)
%    han = plotOverTime(R,dims,type)
%
% Inputs:
%    R - reachSet object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%        for plotting, including added pairs:
%          'Unify', <true/false> compute union of all reachable sets
%          'Set', <whichset> corresponding to
%                   ti ... time-interval reachable set (default)
%                   tp ... time-point reachable set (default if no ti)
%                   y  ... time-interval algebraic set
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, reachSet

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       02-June-2020
% Last update:   15-July-2020 (MW, handling of plot options)
%                01-July-2021 (MP, adding improved unify algorithm)
%                01-April-2023 (MW, speed up fastUnify plotting)
% Last revision: 01-May-2023 (MW, restructure, add unify for tp solutions)
%                12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input
[R,dims,NVpairs,unify,totalsets,whichset] = aux_parseInput(R,varargin{:});

% 2. plot reachable sets
han = aux_plotReachSet(R,dims,NVpairs,unify,totalsets,whichset);

% 3. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [R,dims,NVpairs,unify,totalsets,whichset] = aux_parseInput(R,varargin)
    % default values for the optional input arguments
    dims = setDefaultValues({1},varargin);
    
    % check input arguments
    inputArgsCheck({{R,'att','reachSet','nonempty'};
                    {dims,'att','numeric',{'nonempty','scalar','integer','positive'}}});
    
    % parse input arguments
    NVpairs = readPlotOptions(varargin(2:end),'reachSet');
    [NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',false);
    [NVpairs,totalsets] = readNameValuePair(NVpairs,'UnifyTotalSets','isscalar',1);
    [NVpairs,whichset] = readNameValuePair(NVpairs,'Set','ischar','ti');
    
    % check which set has to be plotted
    whichset = aux_checkSet(R,whichset);
end

function whichset = aux_checkSet(R,whichset)

% must be character vector for switch-expression
if isempty(whichset)
    whichset = '';
end

switch whichset
    case 'ti'
        if isempty(R(1).timeInterval)
            warning("No time-interval reachable set. Time-point reachable set plotted instead.");
            whichset = 'tp';
        end
        
    case 'tp'
        % no issues (should always be computed)

    case 'y'
        if isempty(R(1).timeInterval.algebraic)
            throw(CORAerror('CORA:emptyProperty'));
        end

    otherwise
        % default value
        if isempty(whichset)
            whichset = 'ti';
            if isempty(R(1).timeInterval) 
                if ~isempty(R(1).timePoint)
                    whichset = 'tp';
                else
                    throw(CORAerror('CORA:emptySet'));
                end
            end
        else
            % has to be one of the above three cases...
            throw(CORAerror('CORA:wrongValue','Set',...
                'has to be ''ti'', ''tp'' or ''y''.'));
        end
end

end

function han = aux_plotReachSet(R,dims,NVpairs,unify,totalsets,whichset)
    % plot reachable set

    % check if the reachable sets should be unified to reduce the storage size
    % of the resulting figure
    if unify
        
        if any(strcmp(whichset,{'ti','y'}))
            % check whether fastUnify approach is possible
            if aux_fastUnify(R)
                han = aux_plotFastUnify(R,dims,NVpairs,whichset);
            else
                % unification of polygons
                han = aux_plotUnify(R,dims,NVpairs,whichset);
            end
        elseif strcmp(whichset,'tp')
            % utilize NaN as breaks between lines to plot time-point solutions
            han = aux_plotUnifyTP(R,dims,NVpairs,whichset);
        end
    
    else
        % regular plotting if individual sets
        han = aux_plotStandard(R,dims,NVpairs,whichset);
    
    end
end

function fastUnify = aux_fastUnify(R)
% checks whether fastUnify approach (without polygon union) is possible

% assume true
fastUnify = true;

for i = 1:size(R,1)
    % read out time-interval solution
    Rset = R(i,1).timeInterval;
    
    if ~isempty(Rset)
        % check if intervals are disjoint
        disjoint_check = cellfun(@(x,y){x.supremum-y.infimum},...
            Rset.time(1:end-1),Rset.time(2:end));
        disjoint_check = cell2mat(disjoint_check);
        disjoint_check = unique(disjoint_check);
        if length(disjoint_check) ~= 1 || disjoint_check(1) ~= 0
            fastUnify = false;
            break;
        end
    end
end

end

function [intXmin,intXmax] = aux_getXminmax(set,dims)

if isa(set,'zonotope')
    % faster conversion than calling interval constructor
    % (especially, if there are many sets to plot)
    intXmin = set.c(dims) - sum(abs(set.G(dims,:)));
    intXmax = set.c(dims) + sum(abs(set.G(dims,:)));
else
    intX = interval(project(set,dims));
    intXmax = supremum(intX);
    intXmin = infimum(intX);
end

end

function han = aux_plotFastUnify(R,dims,NVpairs,whichset)

% init lists
x_list = cell(length(R),1);
y_list = cell(length(R),1);

% loop over all reachable sets
for i = 1:size(R,1)
    
    % get desired set
    switch whichset
        case 'ti'
            Rset = R(i,1).timeInterval.set;
            Rtime = R(i,1).timeInterval.time;
        case 'tp'
            Rset = R(i,1).timePoint.set;
            Rtime = R(i,1).timePoint.time;
        case 'y'
            Rset = R(i,1).timeInterval.algebraic;
            Rtime = R(i,1).timeInterval.time;
    end

    x_list{i,1} = zeros(1,4*length(Rset));
    y_list{i,1} = zeros(1,4*length(Rset));
    idx = 1;

    for j = 1:length(Rset)

        % get intervals
        [intXmin,intXmax] = aux_getXminmax(Rset{j},dims);
    
        % use fast unification plotting algorithm
        % add coordinates of interval corners into lists, while
        % iterating through corners of individual intervals
        % clockwise (upper-left corner, upper-right corner,
        % lower-right corner, lower-left corner);
        % add new polygons in the middle to maintain order
        % (up_corner_p1 up_corner_p2 low_corner_p2 low_corner_p1)
        x_list{i}(idx:idx+1) = [infimum(Rtime{j}),supremum(Rtime{j})];
        x_list{i}(end-idx:end-idx+1) = [supremum(Rtime{j}),infimum(Rtime{j})];
        y_list{i}(idx:idx+1) = [intXmax,intXmax];
        y_list{i}(end-idx:end-idx+1) = [intXmin,intXmin];
        % shift index
        idx = idx + 2;
    end
end

% merge lists
x_merged = x_list{1}; y_merged = y_list{1};
for i=2:size(R,1)
    x_merged = [x_merged(1:length(x_merged)/2), x_list{i}, ...
        x_merged(length(x_merged)/2+1:end)];
    y_merged = [y_merged(1:length(y_merged)/2), y_list{i}, ...
        y_merged(length(y_merged)/2+1:end)];
end


% plot patch

if ~ishold
    % patch ignores hold status, thus correcting it manually...
    plot(nan,nan,'HandleVisibility','off')
end

% (list already ordered, skip instantiation of polygon
% which checks for shape simplification)

[NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor',{},CORAcolor("CORA:next"));
han = plotPolygon([x_merged;y_merged],NVpairs{:},'FaceColor',facecolor);

end

function han = aux_plotUnify(R,dims,NVpairs,whichset)

% init polygon
pgon = [];
% turn of warnings (occur in polygon union)
warOrig = warning;
warning('off','all');

% loop over all reachable sets
for i = 1:size(R,1)
    
    % get desired set
    switch whichset
        case 'ti'
            Rset = R(i,1).timeInterval.set;
            Rtime = R(i,1).timeInterval.time;
        case 'tp'
            Rset = R(i,1).timePoint.set;
            Rtime = R(i,1).timePoint.time;
        case 'y'
            Rset = R(i,1).timeInterval.algebraic;
            Rtime = R(i,1).timeInterval.time;
    end

    for j = 1:length(Rset)

        % get intervals
        [intXmin,intXmax] = aux_getXminmax(Rset{j},dims);
    
        % convert to polygon and unite with previous sets
        V = [infimum(Rtime{j}),infimum(Rtime{j}),...
            supremum(Rtime{j}),supremum(Rtime{j});...
            intXmin,intXmax,intXmax,intXmin];
        temp = polygon(V(1,:),V(2,:));
        pgon = pgon | temp;
    end

end

% plot the resulting set
han = plot(pgon,[1,2],NVpairs{:});

% reset original warning status
warning(warOrig);

end

function han = aux_plotUnifyTP(R,dims,NVpairs,whichset)
% unify time-point solutions over time by introducing breaks via NaN
% (note: whichset only as input argument for consistency with other aux_)

% init time and state
t = []; x = [];

% loop over all reachable sets
for i = 1:size(R,1)

    % get desired set
    Rset = R(i,1).timePoint.set;
    Rtime = R(i,1).timePoint.time;

    % create time vector: [t_0, t_0, NaN, t_1, t_1, NaN, ...]
    t_ = cell2mat(Rtime);
    t_ = repmat(t_,1,2);
    t_ = [t_,NaN(length(t_),1)]';
    t_ = t_(:);

    % create state vector: [xmin0, xmax0, NaN, xmin1, xmax1, NaN, ...]
    nrSets = length(Rset);
    x_ = zeros(nrSets*3,1);
    for j = 1:nrSets
        % get bounds for set
        [intXmin,intXmax] = aux_getXminmax(Rset{j},dims);
        x_((j-1)*3+1:(j-1)*3+3) = [intXmin;intXmax;NaN];
    end

    % append to full solution
    t = [t; t_];
    x = [x; x_];
end

% use 'simResult' to convert 'FaceColor'/'EdgeColor' to 'Color' for line
% plot below
NVpairs = readPlotOptions(NVpairs,'simResult');

% plot line
han = plot(t,x,NVpairs{:});

end

function han = aux_plotStandard(R,dims,NVpairs,whichset)
% plot reachable sets as they are

sets = {};

% loop over all reachable sets
for i = 1:size(R,1)
    
    % get desired set
    switch whichset
        case 'ti'
            % fall back to time-point set (hybrid, instant transitions)
            if isempty(R(i,1).timeInterval)
                Rset = R(i,1).timePoint.set;
                Rtime = R(i,1).timePoint.time;
            else
                Rset = R(i,1).timeInterval.set;
                Rtime = R(i,1).timeInterval.time;
            end
        case 'tp'
            Rset = R(i,1).timePoint.set;
            Rtime = R(i,1).timePoint.time;
        case 'y'
            Rset = R(i,1).timeInterval.algebraic;
            Rtime = R(i,1).timeInterval.time;
    end
    
    for j = 1:length(Rset)

        % get intervals
        intX = interval(project(Rset{j},dims));
        if ~isa(Rtime{j},'interval')
            % time-point solution
            intT = interval(Rtime{j});
        else
            intT = Rtime{j};
        end
        I = cartProd(intT,intX);

        % save set
        sets{end+1} = I;
    end
end

han = plotMultipleSetsAsOne(sets,[1,2],NVpairs);

end

% ------------------------------ END OF CODE ------------------------------
