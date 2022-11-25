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
    type = {'b'};
    if nargin >= 2
        dims = varargin{1}; 
    end
    if nargin >= 3
        [linespec,type] = readPlotOptions(varargin(2:end));
        [type,~] = readNameValuePair(type,'Filled','islogical',false);
        type = [{linespec},type];
    end

    % check dimension
    if length(dims) < 2
        error('At least 2 dimensions have to be specified!');
    elseif length(dims) > 3
        error('Only up to 3 dimensions can be plotted!');
    end

    % different types of level sets
    if strcmp(obj.compOp,'==')

        % different methods for the different dimensions
        if length(dims) == 2
            han = plot2Dcontour(obj,dims,type);
        else
            [res,ind] = isSolvable(obj,dims);

            if res
                han = plot3Dsolveable(obj,dims,ind,type);
            else
                han = plot3Dgrid(obj,dims,type); 
            end
        end

    else
        % different methods for the differnt dimensions
        if length(dims) == 2
            han = plot2Dgrid(obj,dims,type);
        else
            han = plot3Dgrid(obj,dims,type); 
        end 
    end
end


% Auxiliary Functions -----------------------------------------------------

function han = plot2Dcontour(obj,dims,type)
% plot 2D level set using Matlabs contour plot function

    % get limits of figure
    ax = gca;
    xlim = get(ax,'Xlim');
    ylim = get(ax,'Ylim');
    
    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);

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
end

function han = plot2Dgrid(obj,dims,type)
% plot 2D level set by gridding the plot area 

    % get limits of figure
    ax = gca;
    xlim = get(ax,'Xlim');
    ylim = get(ax,'Ylim');

    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);
    
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

function han = plot3Dsolveable(obj,dims,ind,type)
% plot 3D level set by solving for one variable

    % get limits of figure
    ax = gca;
    xlim = get(ax,'Xlim');
    ylim = get(ax,'Ylim');
    zlim = get(ax,'Zlim');
    lim = [xlim;ylim;zlim];
    
    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);
    
    % get function handle
    f = obj.solved{ind}.funHan{1}.eq;
    
    % generate grid
    N = 100; 
    ind_ = setdiff(1:3,ind);
    
    dx = (lim(ind_(1),2)-lim(ind_(1),1))/N;
    dy = (lim(ind_(2),2)-lim(ind_(2),1))/N;
    
    dx_ = dx/2; dy_ = dy/2;
    
    x = lim(ind_(1),1)+dx_:dx:lim(ind_(1),2)-dx_;
    y = lim(ind_(2),1)+dy_:dy:lim(ind_(2),2)-dy_;
    
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            p_ = p;
            p_(dims(ind_)) = [X(i,j);Y(i,j)];
            Z(i,j) = f(p_);
        end
    end
    
    if ind == 1
        han = surf(Z,X,Y,'FaceColor',type{1},type{2:end});
    elseif ind == 2
        han = surf(X,Z,Y,'FaceColor',type{1},type{2:end});
    else
        han = surf(X,Y,Z,'FaceColor',type{1},type{2:end});
    end

end

function han = plot3Dgrid(obj,dims,type)
% plot 3D level set by gridding the plot area 

    % get limits of figure
    ax = gca;
    xlim = get(ax,'Xlim');
    ylim = get(ax,'Ylim');
    zlim = get(ax,'Zlim');

    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);
    
    % generate grid
    N = 20;
    dx = (xlim(2)-xlim(1))/N;
    dy = (ylim(2)-ylim(1))/N;
    dz = (zlim(2)-zlim(1))/N;
    
    dx_ = dx/2; dy_ = dy/2; dz_ = dz/2;
    
    x = xlim(1)+dx_:dx:xlim(2)-dx_;
    y = ylim(1)+dy_:dy:ylim(2)-dy_;
    z = zlim(1)+dz_:dz:zlim(2)-dz_;
    
    [X,Y] = meshgrid(x,y);
    
    hold on
    
    % plot all grid cells belonging to the set
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            for k = 1:length(z)
            
                % construct grid cell
                c = p; d = p;
                c(dims) = [X(i,j); Y(i,j); z(k)];
                d(dims) = [dx_; dy_; dz_];
                
                int = interval(c-d,c+d);
                
                % plot the grid cell if it intersects the level set
                if isIntersecting(obj,int,'approx')
                    han = plot(int,dims,type{:},'Filled',true);
                end
            end
        end
    end
end

function [res,ind] = isSolvable(obj,dims)
% check if the level set equation is solvable for one variable

    res = 0;
    
    if ~obj.solveable
       return; 
    end

    % find variable for which the equations are solved
    ind = 0;
    
    for i = 1:length(dims)
        if length(obj.solved{dims(i)}.eq) == 1
            ind = i;
            break;
        end
    end
    
    if ind ~= 0
       res = 1;
    end
end

%------------- END OF CODE --------------