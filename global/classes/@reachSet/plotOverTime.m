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

% Author:       Niklas Kochdumper
% Written:      02-June-2020
% Last update:  15-July-2020 (MW, handling of plot options)
%               01-July-2021 (MP, adding improved unify algorithm)
%               01-April-2023 (MW, speed up fastUnify plotting)
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{R,'att',{'reachSet'},{''}};
                {dims,'att',{'numeric'},{'nonempty','scalar','integer','positive'}}});

% parse input arguments
NVpairs = readPlotOptions(varargin(2:end),'reachSet');
[NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',false);
[NVpairs,whichset] = readNameValuePair(NVpairs,'Set','ischar','ti');

% check which set has to be plotted
whichset = checkSet(R,whichset);

% check if the reachable sets should be unified to reduce the storage size
% of the resulting figure
if unify
    
    pgon = [];
    warOrig = warning;
    warning('off','all');
    
    % flag checking if the fast unified plotting algorithm can be used 
    % (only enabled if the time intervals are disjoint)
    
    fastUnify = false;
        
    if any(strcmp(whichset,{'ti','y'}))
        
        % assume true, check for counterexample
        fastUnify = true;
        
        for i = 1:size(R,1)

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

        % cell-array for each branch of the reachSet object (merged later)
        if fastUnify
            x_list = cell(length(R),1);
            y_list = cell(length(R),1);
        end
    end
    
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

        % init array for data points (4 per set) for i-th branch
        if fastUnify
            x_list{i,1} = zeros(1,4*length(Rset));
            y_list{i,1} = zeros(1,4*length(Rset));
            idx = 1;
        end
        
        for j = 1:length(Rset)

            % get intervals
            if isa(Rset{j},'zonotope')
                % faster conversion than calling interval constructor
                % (especially, if there are many sets to plot)
                intXmin = Rset{j}.Z(dims,1) - sum(abs(Rset{j}.Z(dims,2:end)));
                intXmax = Rset{j}.Z(dims,1) + sum(abs(Rset{j}.Z(dims,2:end)));
            else
                intX = interval(project(Rset{j},dims));
            end
            
            % check flag
            if fastUnify
                % use fast unification plotting algorithm
                % add coordinates of interval corners into lists, while
                % iterating through corners of individual intervals
                % clockwise (upper-left corner, upper-right corner,
                % lower-right corner, lower-left corner);
                % add new polygons in the middle to maintain order
                % (up_corner_p1 up_corner_p2 low_corner_p2 low_corner_p1)
                x_list{i}(idx:idx+1) = [infimum(Rtime{j}),supremum(Rtime{j})];
                x_list{i}(end-idx:end-idx+1) = [supremum(Rtime{j}),infimum(Rtime{j})];
                if isa(Rset{j},'zonotope')
                    y_list{i}(idx:idx+1) = [intXmax,intXmax];
                    y_list{i}(end-idx:end-idx+1) = [intXmin,intXmin];
                else
                    y_list{i}(idx:idx+1) = [supremum(intX),supremum(intX)];
                    y_list{i}(end-idx:end-idx+1) = [infimum(intX),infimum(intX)];
                end
                % shift index
                idx = idx + 2;
                
            % use regular unification plot algorithm
            else
            	% convert to polygon and unite with previous sets
                if isa(Rset{j},'zonotope')
                    V = [infimum(Rtime{j}),infimum(Rtime{j}),...
                        supremum(Rtime{j}),supremum(Rtime{j});...
                        intXmin,intXmax,intXmax,intXmin];
                else
                    V = [infimum(Rtime{j}),infimum(Rtime{j}),...
                        supremum(Rtime{j}),supremum(Rtime{j});...
                        infimum(intX),supremum(intX),...
                        supremum(intX),infimum(intX)];
                end
                temp = polygon(V(1,:),V(2,:));
                pgon = pgon | temp;
            end

        end  
    end
    
    if fastUnify
        % merge lists
        x_merged = x_list{1}; y_merged = y_list{1};
        for i=2:size(R,1)
            x_merged = [x_merged(1:length(x_merged)/2), x_list{i}, ...
                x_merged(length(x_merged)/2+1:end)];
            y_merged = [y_merged(1:length(y_merged)/2), y_list{i}, ...
                y_merged(length(y_merged)/2+1:end)];
        end
        % plot patch (list already ordered, skip instantiation of polygon
        % which checks for shape simplification)
        % linespec 'b' overwritten by NVpairs, but patch requires it...
        han = patch(x_merged,y_merged,'b',NVpairs{:});
    else
    
        % plot the resulting set
        han = plot(pgon,[1 2],NVpairs{:});
    end
    
    warning(warOrig);
    
else
    hold on;

    % save color index
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;
    
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
            intX = interval(project(Rset{j},dims));
            if ~isa(Rtime{j},'interval')
                % time-point solution
                intT = interval(Rtime{j});
            else
                intT = Rtime{j};
            end
            int = cartProd(intT,intX);

            % plot interval
            han_ij = plot(int,[1,2],NVpairs{:});

            if i == 1 && j == 1
                han = han_ij;
                % don't display subsequent plots in legend
                NVpairs = [NVpairs, {'HandleVisibility','off'}];
            end
        end
    end

    % correct color index
    updateColorIndex(oldColorIndex);
end

if nargout == 0
    clear han;
end

end


% Auxiliary function ------------------------------------------------------

function whichset = checkSet(R,whichset)

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

%------------- END OF CODE --------------