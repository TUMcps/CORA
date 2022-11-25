function han = plot(R,varargin)
% plot - plots a 2-dimensional projection of the reachable set
%
% Syntax:  
%    han = plot(R)
%    han = plot(R,dim)
%    han = plot(R,dim,plotOptions)
%
% Inputs:
%    R - reachSet object
%    dims - (optional) dimensions that should be projected
%    plotOptions - (optional) LineSpecification and Name-Value pairs
%        for plotting, including added pairs:
%          'Filled' <filled plot or not>
%          'Order', <zonotope order>,
%          'Splits', <number of splits for refinement of polyZonotopes>
%          'Unify', <compute union of all reachable sets>
%
% Outputs:
%    han - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-June-2020
% Last update:  15-July-2020 (MW, merge with plotFilled, plotOptions)
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dims = [1,2];
linespec = 'b';
NVpairs = {};
order = [];
splits = [];
unify = 0;

% parse input arguments
if nargin > 1 && ~isempty(varargin{1})
   dims = varargin{1}; 
end
if nargin > 2 && ~isempty(varargin{2})
    plotOptions = varargin(2:end);
    [linespec,NVpairs] = readPlotOptions(plotOptions);
    [NVpairs,order] = readNameValuePair(NVpairs,'Order','isscalar');
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits','isscalar');
    [NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',unify);
    
    % in reachSet name-value pair 'Filled',true|false is default true,
    % this must be added because called plot functions are default false
    [NVpairs,filled] = readNameValuePair(NVpairs,'Filled','islogical');
    if isempty(filled) || (islogical(filled) && filled)
        NVpairs = [NVpairs,'Filled',true];
    end
end
% check dimension
if length(dims) < 2
    error('At least 2 dimensions have to be specified!');
elseif length(dims) > 3
    error('Only up to 3 dimensions can be plotted!');
elseif length(dims) == 3
    unify = 0;
end


% check if the reachable sets should be unified to reduce the storage size
% of the resulting figure
if unify

    pgon = [];
    warOrig = warning;
    warning('off','all');
    
    % loop over all reachable sets
    for i = 1:size(R,1)
        
        % get reachable set
        if ~isempty(R(i,1).timeInterval)
            set = R(i,1).timeInterval.set;
        else
            set = R(i,1).timePoint.set;
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
                   isa(temp,'mptPolytope') || isa(temp,'zonoBundle') || ...
                   isa(temp,'conZonotope')
                V = vertices(temp);
                temp = polygon(V(1,:),V(2,:));
            else
                error('Setting "Unify" is not supported for this set representation!');
            end
            
            % unite all polygons
            pgon = pgon | temp;
        end
    end
    
    % plot the resulting set
    han = plot(pgon,[1 2],linespec,NVpairs{:});
    
    warning(warOrig);
    
else
    
    % loop over all reachable sets
    hold on

    for i = 1:size(R,1)

        % get reachable set
        if ~isempty(R(i,1).timeInterval)
            set = R(i,1).timeInterval.set;
        else
            set = R(i,1).timePoint.set;
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
                    han = plot(zonotope(temp),dims_,linespec,NVpairs{:});
                else
                    han = plot(temp,dims_,linespec,'Splits',splits,NVpairs{:});
                end
            else
                han = plot(temp,dims_,linespec,NVpairs{:});
            end
        end
    end
end

%------------- END OF CODE --------------