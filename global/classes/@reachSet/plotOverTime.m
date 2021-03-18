function han = plotOverTime(R,varargin)
% plotOverTime - plots the reachable set over time
%
% Syntax:  
%    han = plotOverTime(R)
%    han = plotOverTime(R,dim,linespec)
%    han = plotOverTime(R,dim,linespec,NVpairs)
%
% Inputs:
%    R - reachSet object
%    dim - dimension that should be projected
%    linespec - Line Specification properties
%    NVpairs - name-value pairs (Patch Properties)
%
% Outputs:
%    han - handle for the resulting graphic object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, reachSet

% Author:       Niklas Kochdumper
% Written:      02-June-2020
% Last update:  15-July-2020 (MW, handling of plot options)
% Last revision:---

%------------- BEGIN CODE --------------

% default values for the optional input arguments
dim = 1;
linespec = 'b';
unify = 0;
NVpairs = {};

% parse input arguments
if nargin > 1 && ~isempty(varargin{1})
    dim = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    plotOptions = varargin(2:end);
    [linespec,NVpairs] = readPlotOptions(plotOptions);
    [NVpairs,unify] = readNameValuePair(NVpairs,'Unify','islogical',unify);
end

% check if the reachable sets should be unified to reduce the storage size
% of the resulting figure
if unify
    
    pgon = [];
    warOrig = warning;
    warning('off','all');
    
    % loop over all reachable sets
    for i = 1:size(R,1)
        
        if ~isempty(R(i,1).timeInterval)
            Rset = R(i,1).timeInterval.set;
            Rtime = R(i,1).timeInterval.time;
        else
            Rset = R(i,1).timePoint.set;
            Rtime = R(i,1).timePoint.time;
        end
        
        for j = 1:length(Rset)

            % get intervals
            intX = interval(project(Rset{j},dim));
            intT = Rtime{j};

            int = cartProd(intT,intX);

            % convert to polygon and unite with previous sets
            V = [infimum(int(1)),infimum(int(1)),supremum(int(1)),supremum(int(1)); ...
                infimum(int(2)),supremum(int(2)),supremum(int(2)),infimum(int(2))];
            temp = polygon(V(1,:),V(2,:));
            pgon = pgon | temp;
        end
    end
    
    % plot the resulting set
    han = plot(pgon,[1 2],linespec,NVpairs{:},'Filled',true);
    
    warning(warOrig);
    
else
    hold on;
    
    % loop over all reachable sets
    for i = 1:size(R,1)
        
        if ~isempty(R(i,1).timeInterval)
            Rset = R(i,1).timeInterval.set;
            Rtime = R(i,1).timeInterval.time;
        else
            Rset = R(i,1).timePoint.set;
            Rtime = R(i,1).timePoint.time;
        end
        
        for j = 1:length(Rset)

            % get intervals
            intX = interval(project(Rset{j},dim));
            intT = Rtime{j};

            int = cartProd(intT,intX);

            % plot interval
            han = plot(int,[1,2],linespec,NVpairs{:},'Filled',true);
        end
    end
end

%------------- END OF CODE --------------