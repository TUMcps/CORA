function han = plot(R,varargin)
% plot - plots a projection of the reachable set
%
% Syntax:  
%    han = plot(R)
%    han = plot(R,dims)
%    han = plot(R,dims,plotOptions)
%
% Inputs:
%    R - reachSet object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs),
%        including added pairs:
%          'Order', <order> (polynomial) zonotope order
%          'Splits', <splits> number of splits for refinement
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
% See also: reachSet

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-June-2020
% Last update:  15-July-2020 (MW, merge with plotFilled, plotOptions)
%               29-October-2021 (MW, remove linespec, 'Filled')
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = setDefaultValues({[1,2]},varargin);

% check input arguments
inputArgsCheck({{R,'att',{'reachSet'},{''}};
                {dims,'att',{'numeric'},{'nonempty','vector','integer','positive'}}});

% check name-value pairs
NVpairs = readPlotOptions(varargin(2:end),'reachSet');
[NVpairs,order] = readNameValuePair(NVpairs,'Order','isscalar');
[NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
[NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',false);
[NVpairs,whichset] = readNameValuePair(NVpairs,'Set','ischar','ti');

% check dimension
if length(dims) < 2
    throw(CORAerror('CORA:plotProperties',2));
elseif length(dims) > 3
    throw(CORAerror('CORA:plotProperties',3));
elseif length(dims) == 3
    unify = false;
end

% check which set has to be plotted
whichset = checkSet(R,whichset);

% check if the reachable sets should be unified to reduce the storage size
% of the resulting figure
if unify

    pgon = [];
    warOrig = warning;
    warning('off','all');
    
    % loop over all reachable sets
    for i = 1:size(R,1)
        
        % get desired set
        switch whichset
            case 'ti'
                set = R(i,1).timeInterval.set;
            case 'tp'
                set = R(i,1).timePoint.set;
            case 'y'
                set = R(i,1).timeInterval.algebraic;
        end
        
        % loop over all time steps
        for j = 1:length(set)

            % project set to desired dimensions
            temp = project(set{j},dims);

            % order reduction
            if ~isempty(order)
                temp = reduce(temp,'girard',order);
            end
           
            % convert set to polygon object
            if isa(temp,'polyZonotope') || isa(temp,'conPolyZono')
                if isempty(splits)
                    V = vertices(zonotope(temp));
                    temp = polygon(V(1,:),V(2,:));
                else
                    temp = polygon(temp,splits);   
                end
            elseif isa(temp,'zonotope') || isa(temp,'interval') || ...
                   isa(temp,'mptPolytope') || isa(temp,'conZonotope')
                V = vertices(temp);
                temp = polygon(V(1,:),V(2,:));
            elseif isa(temp,'zonoBundle')                
                % compute polytopes
                P = cell(temp.parallelSets,1);
                for p = 1:temp.parallelSets
                    % delete zero-length generators
                    Z = deleteZeros(temp.Z{p});
                    % convert to polyshape (Matlab built-in class)
                    temp_poly = polygon(Z);
                    V = temp_poly(:,2:end);
                    P{p} = polyshape(V(1,:),V(2,:));
                end

                % intersect polytopes
                Pint = P{1};
                for p = 2:temp.parallelSets
                    Pint = intersect(Pint,P{p});
                end

                % get vertices, init polygon
                V = Pint.Vertices';
                temp = polygon(V(1,:),V(2,:));
            else
                throw(CORAerror('CORA:specialError',['Setting "Unify" is '...
                    'not supported for this set representation!']));
            end
            
            % unite all polygons
            pgon = pgon | temp;
        end
    end
    
    % plot the resulting set
    han = plot(pgon,[1 2],NVpairs{:});
    
    warning(warOrig);
    
else
    
    % loop over all reachable sets
    hold on

    % save color index
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;

    for i = 1:size(R,1)

        % get desired set
        switch whichset
            case 'ti'
                set = R(i,1).timeInterval.set;
            case 'tp'
                set = R(i,1).timePoint.set;
            case 'y'
                set = R(i,1).timeInterval.algebraic;
        end

        % loop over all time steps
        for j = 1:length(set)

            % project set to desired dimensions
            temp = project(set{j},dims);
            
            if length(dims) == 2
               dims_ = [1,2]; 
            else
               dims_ = [1,2,3]; 
            end

            % order reduction
            if ~isempty(order)
                temp = reduce(temp,'girard',order);
            end

            % plot the set
            if isa(temp,'polyZonotope')
                if isempty(splits)
                    han_ij = plot(zonotope(temp),dims_, NVpairs{:});
                else
                    han_ij = plot(temp,dims_,'Splits',splits, NVpairs{:});
                end
            else
                han_ij = plot(temp,dims_,NVpairs{:});
            end

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


% Auxiliary function
function whichset = checkSet(R,whichset)

% must be character vector for switch-expression to work properly
if isempty(whichset)
    % default value
    if ~isempty(R(1).timeInterval)
        whichset = 'ti';
    else
        whichset = 'tp';
    end
end

switch whichset
    case 'ti'
        if isempty(R(1).timeInterval) && ~isempty(R(1).timePoint)
            warning("No time-interval reachable set. Time-point reachable set plotted instead.");
            whichset = 'tp';
        end
        
    case 'tp'
        % no issues (should always be there)

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