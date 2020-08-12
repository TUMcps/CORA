function han = plot(obj,varargin)
% plot - plots 2-dimensional projection of a levelSet
%
% Syntax:  
%    plot(obj)
%    plot(obj,dims)
%    plot(obj,dims,type)
%
% Inputs:
%    obj - levelSet object
%    dims - (optional) dimensions that should be projected (optional).
%               assumption: other entries of the normal vector are zeros
%    type - (optional) plot settings (LineSpec and name-value pairs)
%
% Outputs:
%    han - handle to graphics object
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    figure; hold on; box on;
%    xlim([-3,3]);
%    ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      19-July-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
dims = [1,2];
type = 'b';
if nargin >= 2
    dims = varargin{1}; 
end
if nargin >= 3
	type = varargin(2:end);
end

% get limits of figure
ax = gca;
xlim = get(ax,'Xlim');
ylim = get(ax,'Ylim');

% substitute all remaining entries with zero
p = zeros(obj.dim,1);

% different types of level sets
if strcmp(obj.compOp,'==')
    
    % generate contour plot
    N = 30;
    x = xlim(1):(xlim(2)-xlim(1))/N:xlim(2);
    y = ylim(1):(ylim(2)-ylim(1))/N:ylim(2);

    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));

    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            p_ = p;
            p_(dims) = [X(i,j);Y(i,j)];
            Z(i,j) = obj.funHan(p_);
        end
    end

    [~,han] = contour(X,Y,Z,[0,0],type{:});

else
    
    % generate grid
    N = 200;
    dx = (xlim(2)-xlim(1))/N;
    dy = (ylim(2)-ylim(1))/N;
    
    dx_ = dx/2;
    dy_ = dy/2;
    
    x = xlim(1)+dx_:dx:xlim(2)-dx_;
    y = ylim(1)+dy_:dy:ylim(2)-dy_;
    
    [X,Y] = meshgrid(x,y);
    
    hold on
    
    % plot all grid cells belonging to the set
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            
            % evaluate level set function for the center of the grid cell
            p_ = p;
            p_(dims) = [X(i,j);Y(i,j)];
            val = obj.funHan(p_);
            
            % plot the grid cell
            if val <= 0
                x = X(i,j) + [-dx_ -dx_ dx_ dx_];
                y = Y(i,j) + [-dy_ dy_ dy_ -dy_];
                
                han = fill(x,y,type{:},'EdgeColor','none');
            end
        end
    end  
end

%------------- END OF CODE --------------