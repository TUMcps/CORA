function han = plot(ls,varargin)
% plot - plots a projection of a level set
%
% Syntax:
%    han = plot(ls)
%    han = plot(ls,dims)
%    han = plot(ls,dims,varargin)
%
% Inputs:
%    ls - levelSet object
%    dims - (optional) dimensions for projection
%           (assumption: other entries of the normal vector are zeros)
%    varargin - (optional) plot settings (LineSpec and Name-Value pairs)
%               - 'Splits': number of splits to plot sets with '<=','<'
%               - 'PlotMethod': one of {'outer','inner'}, default: 'outer'
%
% Outputs:
%    han - handle to the graphics object
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    figure; hold on; box on; xlim([-3,3]); ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       19-July-2019
% Last update:   22-May-2023 (TL, speed up plotting of '<=' levelSets)
%                26-July-2023 (TL, getUnboundedAxisLimits)
% Last revision: 12-July-2023 (TL, restructure)

% ------------------------------ BEGIN CODE -------------------------------

% 1. parse input arguments
[ls,dims,NVpairs,splits,plotMethod] = aux_parseInput(ls,varargin{:});

% 2. preprocess
[ls,dims] = aux_preprocess(ls,dims);

% 3. plot n-dimensional set
han = aux_plotNd(ls,dims,NVpairs,splits,plotMethod);

% 4. clear han
if nargout == 0
    clear han;
end

end


% Auxiliary functions -----------------------------------------------------

function [ls,dims,NVpairs,splits,plotMethod] = aux_parseInput(ls,varargin)
    % parse input arguments
    dims = setDefaultValues({[1,2]},varargin);

    % check input arguments
    inputArgsCheck({{ls,'att','levelSet'};
                    {dims,'att','numeric',{'nonnan','integer','vector','positive'}}});

    % check dimension
    if length(dims) < 1
        throw(CORAerror('CORA:plotProperties',1));
    elseif length(dims) > 3
        throw(CORAerror('CORA:plotProperties',3));
    end

    % read out plot options and additional name-value pairs
    NVpairs = varargin(2:end);
    % 'Splits' given?
    [NVpairs,splits] = readNameValuePair(NVpairs,'Splits');
    if isempty(splits)
        if length(dims) == 3
            splits = 4;
        else
            splits = 7;
        end
    end
    % 'PlotMethod' given?
    [NVpairs,plotMethod] = readNameValuePair(NVpairs,'PlotMethod',{},'outer');

    % check name-value pairs
    if CHECKS_ENABLED
        if ~isnumeric(splits) || ~isscalar(splits)
            throw(CORAerror("CORA:wrongValue","name-value pair 'Splits'",'numeric scalar'))
        end
        if ~ismember(plotMethod,{'inner','outer'})
            throw(CORAerror("CORA:wrongValue","name-value pair 'PlotMethod'","'outer','inner'"))
        end
    end
end

function [ls,dims] = aux_preprocess(ls,dims)
    % preprocess

    if length(dims) == 1
        % add artificial dimension at 2nd dimension
        dim_old = 1:length(ls.vars);
        dim_old = 2 + dim_old; % shift
        dim_old(dims) = 1;
        ls = lift_(ls, length(dim_old)+2, dim_old);
        dims = [1;2];
    end
end

function han = aux_plotNd(ls,dims,NVpairs,splits,plotMethod)

    hold on; han = [];

    % different types of level sets
    if all(strcmp(ls.compOp,'=='))      % equality constraints only
    
        % different methods for the different dimensions
        if length(dims) == 2

            if length(ls.eq) == 1
                han = aux_plot2Dcontour(ls,dims,[],NVpairs);
            else
                han = aux_plot2Dsplit(ls,dims,splits,plotMethod,NVpairs);
            end

        elseif length(ls.eq) == 1

            % split equation into product of terms
            lsFac = aux_factorLevelSet(ls);

            for i = 1:length(lsFac)

                [res,ind] = aux_isSolvable(lsFac{i},dims);
                
                if res
                    han = aux_plot3Dsolvable(lsFac{i},dims,ind,NVpairs);
                else
                    han = aux_plot3Dcontour(lsFac{i},dims,[],NVpairs); 
                end
            end

        elseif length(ls.eq) == 2
            han = aux_plot3Dintersection(ls,dims,[],NVpairs);
        else
            han = aux_plot3Dsplit(ls,dims,splits,plotMethod,NVpairs);
        end
    
    elseif all(~strcmp(ls.compOp,'=='))    % inequality constraints only

        % different methods for the different dimensions
        if length(dims) == 2
            han = aux_plot2Dinequality(ls,dims,NVpairs);
        else
            han = aux_plot3Dsplit(ls,dims,splits,plotMethod,NVpairs); 
        end 

    else    % mixed equality and inequality constraints

        % split into equality and inequality constraints
        [lsEq,lsIneq] = splitEquality(ls);
        ineq.set = lsIneq;

        % different methods for different dimensions
        if length(dims) == 2
            if length(lsEq.eq) == 1
                han = aux_plot2Dcontour(lsEq,dims,ineq,NVpairs);
            else
                han = aux_plot2Dsplit(ls,dims,splits,plotMethod,NVpairs);
            end
        else
            if length(lsEq.eq) == 1
                han = aux_plot3Dcontour(lsEq,dims,ineq,NVpairs);
            elseif length(lsEq.eq) == 2
                han = aux_plot3Dintersection(lsEq,dims,ineq,NVpairs);
            else
                han = aux_plot3Dsplit(ls,dims,splits,plotMethod,NVpairs);
            end
        end

    end
end

function han = aux_plot2Dcontour(obj,dims,ineq,NVpairs)
% plot 2D level set using Matlabs contour plot function

    % read plotOptions for contour
    NVpairs = readPlotOptions(NVpairs,'contour');

    % get limits of current plot
    [xLim,yLim,~] = getUnboundedAxisLimits();
    
    % substitute all remaining entries with zero
    p = zeros(obj.dim,1);

    % generate contour plot
    N = 30;
    x = xLim(1):(xLim(2)-xLim(1))/N:xLim(2);
    y = yLim(1):(yLim(2)-yLim(1))/N:yLim(2);

    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));

    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            p_ = p;
            p_(dims) = [X(i,j);Y(i,j)];
            Z(i,j) = obj.funHan(p_);
        end
    end

    % check if additional inequality constraints are given
    if isempty(ineq)

        % level at which contour is plotted: always at z = 0
        level = [0 0];
        [~,han] = contour(X,Y,Z,level,NVpairs{:});
        updateColorIndex; % does not get updated for contour plots

    else

        % compute contour line
        M = contourc(x,y,Z,[0,0]); 

        % split into connected segments
        cont = {};

        while ~isempty(M)
            len = M(2,1); cont{end+1} = M(:,2:1+len); M = M(:,2+len:end);
        end

        % remove points that violate the inequality constraints
        seg = aux_removeSegmentsIneqCons(cont,ineq.set);

        % fix read plot options as plot() is used instead of contour()
        [NVpairs,edgecolor] = readNameValuePair(NVpairs,'edgecolor');
        NVpairs = [NVpairs,{'Color',edgecolor}];

        % plot contour line individually
        for i = 1:length(seg)
            han = plot(seg{i}(1,:),seg{i}(2,:),NVpairs{:});
        end
    end
end

function han = aux_plot2Dinequality(obj,dims,NVpairs)
% plot 2D level set with inequality constraints by filling the correct side
% of the contour lines

    pgon = []; pgon_ = []; han = [];

    % select the correct dimensions
    obj = aux_selectCorrectDimensions(obj,dims);

    % read plotOptions, since always surf ~ fill called
    NVpairs = readPlotOptions(NVpairs,'fill');

    % get limits of current plot
    [xLim,yLim] = getUnboundedAxisLimits();
    lim = [xLim;yLim]; lim = lim(dims,:);

    % setup meshgrid for x-y-plane
    N = 50; 
    
    dx = (lim(1,2)-lim(1,1))/N;
    dy = (lim(2,2)-lim(2,1))/N;
    
    x = lim(1,1):dx:lim(1,2);
    y = lim(2,1):dy:lim(2,2);
    
    [X,Y] = meshgrid(x,y);

    % compute z-coordinates for all grid values
    Z = repmat({zeros(size(X))},[length(obj.eq),1]);
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            p = obj.funHan([X(i,j);Y(i,j)]);
            for k = 1:length(p)
                Z{k}(i,j) = p(k);
            end
        end
    end

    % loop over all equations
    for i = 1:length(obj.eq)

        % construct level set for this equation
        ls = levelSet(obj.eq(i),obj.vars,'==');

        % compute contour line
        M = contourc(x,y,Z{i},[0,0]); 

        % check if the contour is empty
        if isempty(M)
            [pgon_,~] = aux_getPolygons([],ls,[],lim);
        end

        % split into connected segments
        cont = {};

        while ~isempty(M)
            len = M(2,1); cont{end+1} = M(:,2:1+len); M = M(:,2+len:end);
        end

        % compute polygons for the inside region for each segment
        for j = 1:length(cont)

            % compute polygons
            [pgonNew,~] = aux_getPolygons(cont{j},ls,[],lim);

            % combine with previous polygons
            if j == 1
                pgon_ = pgonNew;
            else
                if isIntersecting(pgon_,pgonNew)
                    pgon_ = pgon_ & pgonNew;
                else
                    pgon_ = pgon_ | pgonNew;
                end
            end
        end

        % terminate if polygon is empty
        if representsa(pgon_,'emptySet')
            break;
        end

        % combine with polygons from previous equations
        if i == 1
            pgon = pgon_;
        else
            pgon = pgon & pgon_;
        end
    end

    % plot the resulting polygon
    if ~representsa(pgon,'emptySet')
        han = plot(pgon,dims,NVpairs{:});
    end
end

function han = aux_plot2Dsplit(obj,dims,splits,plotMethod,NVpairs)
% plot 2D level set by splitting the plot area 

    % re-read plotOptions, since always fill called
    NVpairs = readPlotOptions(NVpairs,'fill');


    % get limits of current plot
    [xLim,yLim,~] = getUnboundedAxisLimits();
    space = interval([xLim(1);yLim(1)],[xLim(2);yLim(2)]);

    % determine subspaces to be plotted
    subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod);

    % plot ---

    % read plot settings (TODO: move to new function (also for reachset etc.)?)
    holdStatus = ishold;
    if ~holdStatus
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        set(gca(),'ColorOrderIndex',1);

        % reset limits
        xlim(xLim);
        ylim(yLim);
    end
    hold on;
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;

    % plot sets (TODO: unify)
    for i=1:length(subSpaces)
        han_i = plot(subSpaces{i},[1,2],NVpairs{:});
        if i == 1
            han = han_i;
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
    if isempty(subSpaces)
        [NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');
        han = fill(nan,nan,facecolor,NVpairs{:});
    end

    % restore plot settings
    updateColorIndex(oldColorIndex);
    if ~holdStatus
        hold off;
    end
end

function subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod)
    % refine space and return included interval subspaces
    subSpaces = {};

    % init n-dim space (with all other dimenions = 0)
    nSpace = interval(zeros(dim(obj),1));
    nSpace(dims) = space;

    % range bounding using fast interval arithmetic
    val = obj.funHan(nSpace);

    if all(val.sup <= 0) % compOp always '<=' or '<' here
        % if all equations satisfied
        subSpaces{1} = space;
    elseif any(val.inf > 0)
        % if any equation is unsatisfiable
        % don't include subspace
    else
        % not (yet) decidedable
        if splits == 0
            switch plotMethod
                case 'outer'
                    % plot
                    subSpaces{1} = space;
                case 'inner'
                    % don't plot
            end
        else
            % refine subspace
            partSpace = partition(space, 2);
            
            % check subspaces
            for i = 1:length(partSpace)
                subSpaces_i = aux_refineSpace(obj,dims,partSpace{i},splits-1,plotMethod);
                subSpaces = [subSpaces,subSpaces_i];
            end
        end
    end
end

function han = aux_plot3Dsolvable(obj,dims,ind,NVpairs)
% plot 3D level set by solving for one variable

    % re-read plotOptions, since always surf ~ fill called
    NVpairs = readPlotOptions(NVpairs,'fill');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    lim = [xLim;yLim;zLim];
    
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
    
    % setup meshgrid
    [X,Y] = meshgrid(x,y);
    Z = zeros(size(X));
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            p_ = p;
            p_(dims(ind_)) = [X(i,j);Y(i,j)];
            Z(i,j) = f(p_);
        end
    end

    % replace numbers outside of the range with NaN to prevent them from
    % getting plottet
    Z_ = inf*ones(size(Z,1)+2,size(Z,2)+2);
    Z_(2:end-1,2:end-1) = Z;

    for i = 2:size(Z_,1)-1
        for j = 2:size(Z_,2)-1
            tmp = Z_([i-1,i,i+1],[j-1,j,j+1]);
            if all(all(tmp < lim(ind,1) | tmp > lim(ind,2)))
                Z(i-1,j-1) = NaN;
            end
        end
    end
    
    % plot
    if ind == 1
        han = surf(Z,X,Y,NVpairs{:});
    elseif ind == 2
        han = surf(X,Z,Y,NVpairs{:});
    else
        han = surf(X,Y,Z,NVpairs{:});
    end

end

function han = aux_plot3Dsplit(obj,dims,splits,plotMethod,NVpairs)
% plot 3D level set by splitting the plot area

    % re-read plotOptions, since always fill called
    NVpairs = readPlotOptions(NVpairs,'fill');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    space = interval([xLim(1);yLim(1);zLim(1)],[xLim(2);yLim(2);zLim(2)]);

    % determine subspaces to be plotted
    subSpaces = aux_refineSpace(obj,dims,space,splits,plotMethod);

    % plot ---

    % read plot settings (TODO: move to new function (also for reachset etc.)?)
    holdStatus = ishold;
    if ~holdStatus
        plot(NaN,NaN,'HandleVisibility','off');
        % reset color index (before readPlotOptions!)
        set(gca(),'ColorOrderIndex',1);

        % reset limits
        xlim(xLim);
        ylim(yLim);
        zlim(zLim);
    end
    hold on;
    ax = gca();
    oldColorIndex = ax.ColorOrderIndex;

    % plot sets
    for i=1:length(subSpaces)
        han_i = plot(subSpaces{i},[1,2,3],NVpairs{:});
        if i == 1
            han = han_i;
            NVpairs = [NVpairs, {'HandleVisibility','off'}];
        end
    end
    if isempty(subSpaces)
        [NVpairs,facecolor] = readNameValuePair(NVpairs,'FaceColor');
        han = fill3(nan,nan,nan,facecolor,NVpairs{:});
    end

    % restore plot settings
    updateColorIndex(oldColorIndex);
    if ~holdStatus
        hold off;
    end
end

function [res,ind] = aux_isSolvable(obj,dims)
% check if the level set equation is solvable for one variable

    res = false; ind = 0;
    
    if ~obj.solvable
       return; 
    end

    % find variable for which the equations are solved    
    for i = 1:length(dims)
        if length(obj.solved{dims(i)}.eq) == 1
            ind = i;
            break;
        end
    end
    
    if ind ~= 0
       res = true;
    end
end

function han = aux_plot3Dcontour(obj,dims,ineq,NVpairs)
% plot 3D level set by interpolating between the contour lines at different
% z-levels

    han = [];

    % select the correct dimensions
    obj = aux_selectCorrectDimensions(obj,dims);

    % re-read plotOptions, since always surf ~ fill called
    NVpairs = readPlotOptions(NVpairs,'fill');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    lim = [xLim;yLim;zLim]; lim = lim(dims,:);

    % compute triangular faces that represent the level set
    F = aux_getFaces3D(obj,lim);

    % remove parts of the surfaces that violate the inequality constraints
    if ~isempty(ineq)
        for j = 1:length(F)

            F_ = {};

            for k = 1:length(F{j})
                tmp = contains(ineq.set,F{j}{k});
                if all(tmp)
                    F_{end+1,1} = F{j}{k};
                elseif any(tmp)
                    tmp = aux_intersect3Dinequality(F{j}{k},ineq.set);
                    F_ = [F_;tmp];
                end
            end

            F{j} = F_;
        end
    end

    % plot the surfaces
    for j = 1:length(F)
        for k = 1:length(F{j})
            han = fill3(F{j}{k}(1,:),F{j}{k}(2,:),F{j}{k}(3,:),'r',NVpairs{:});
        end
    end
end

function han = aux_plot3Dintersection(obj,dims,ineq,NVpairs)
% plot intersection of two 3D level sets

    han = [];

    % select the correct dimensions
    obj = aux_selectCorrectDimensions(obj,dims);

    % re-read plotOptions, since always surf ~ fill called
    NVpairs = readPlotOptions(NVpairs,'surf');

    % get limits of current plot
    [xLim,yLim,zLim] = getUnboundedAxisLimits();
    lim = [xLim;yLim;zLim]; lim = lim(dims,:);
    dims = 1:3;

    % detect case where level sets are identical
    ls1 = levelSet(obj.eq(1),obj.vars,"==");
    lsFac1 = aux_factorLevelSet(ls1); 
    ind1 = 1:length(lsFac1);

    ls2 = levelSet(obj.eq(2),obj.vars,"==");
    lsFac2 = aux_factorLevelSet(ls2);
    ind2 = 1:length(lsFac2);

    lsEq = {};

    for i = 1:length(lsFac1)
        for j = 1:length(lsFac2)
            if isequal(lsFac1{i}.eq,lsFac2{j}.eq)
                lsEq{end+1} = lsFac1{i};
                ind1 = setdiff(ind1,i);
                ind2 = setdiff(ind2,j);
            end
        end
    end

    lsFac1 = lsFac1(ind1); lsFac2 = lsFac2(ind2);

    % compute triangular faces that represent the first level set
    F1 = [];

    for i = 1:length(lsFac1)
        
        [res,ind] = aux_isSolvable(lsFac1{i},dims);

        if res && lsFac1{i}.solvable && lsFac1{i}.solved{3}.solvable && ...
                        length(lsFac1{i}.solved{3}.eq) == 1 && ...
                        isempty(symvar(lsFac1{i}.solved{3}.cond))
            ind = 3;
        end

        if res
            F1_ = aux_getFaces3Dsolvable(lsFac1{i},lim,ind);
        else
            F1_ = aux_getFaces3D(lsFac1{i},lim);
        end

        F1 = aux_combineFaces(F1,F1_);
    end

    % compute triangular faces that represent the second level set
    F2 = [];

    for i = 1:length(lsFac2)
        
        [res,ind] = aux_isSolvable(lsFac2{i},dims);

        if res && lsFac2{i}.solvable && lsFac2{i}.solved{3}.solvable && ...
                        length(lsFac2{i}.solved{3}.eq) == 1 && ...
                        isempty(symvar(lsFac2{i}.solved{3}.cond))
            ind = 3;
        end

        if res
            F2_ = aux_getFaces3Dsolvable(lsFac2{i},lim,ind);
        else
            F2_ = aux_getFaces3D(lsFac2{i},lim);
        end

        F2 = aux_combineFaces(F2,F2_);
    end

    % compute intersection of the faces for the two level sets
    lines = {}; faces = {};

    for i = 1:length(F1)
        for j = 1:length(F1{i})

            if ~isempty(F2)
                for k = 1:length(F2{i})
    
                    [li,fa] = aux_intersectTriangularSurfaces(F1{i}{j}, ...
                                                                F2{i}{k});
    
                    if ~isempty(li)
                        lines{end+1} = li;
                    end
    
                    if ~isempty(fa)
                        faces{end+1} = fa;
                    end
                end
            end
        end
    end

    % remove line segments that violate the inequality constraints
    if ~isempty(ineq)
        lines = aux_removeSegmentsIneqCons(lines,ineq.set);
    end

    % fix read plot options as plot3() is used instead of surf()
    [NVpairs,edgecolor] = readNameValuePair(NVpairs,'edgecolor');
    NVpairs = [NVpairs,{'Color',edgecolor}];

    % plot the lines for the intersections
    for j = 1:length(lines)
        han = plot3(lines{j}(1,:),lines{j}(2,:),lines{j}(3,:),'r',NVpairs{:});
    end

    % plot common constraints between the two level sets
    [NVpairs,~] = readNameValuePair(NVpairs,'facecolor');
    NVpairs = [NVpairs,{'FaceColor',edgecolor}];

    for i = 1:length(lsEq)
        [res,ind] = aux_isSolvable(lsEq{i},dims);
                    
        if res && isempty(ineq)
            han = aux_plot3Dsolvable(lsEq{i},dims,ind,NVpairs);
        else
            han = aux_plot3Dcontour(lsEq{i},dims,ineq,NVpairs); 
        end
    end

    % plot faces
    NVpairs = readPlotOptions(NVpairs,'fill');

    for i = 1:length(faces)
        han = fill3(faces{i}(1,:),faces{i}(2,:),faces{i}(3,:),'r',NVpairs{:});
    end
end

function F = aux_getFaces3D(ls,lim)
% compute the faces that represent an equality constraint in 3D

    % setup meshgrid for x-y-plane
    N = 50; 
    
    dx = (lim(1,2)-lim(1,1))/N;
    dy = (lim(2,2)-lim(2,1))/N;
    dz = (lim(3,2)-lim(3,1))/N;
    
    x = lim(1,1):dx:lim(1,2);
    y = lim(2,1):dy:lim(2,2);
    z = lim(3,1):dz:lim(3,2);
    
    [X,Y] = meshgrid(x,y);

    % compute z-coordinates for all grid values
    Z = repmat({zeros(size(X))},[length(z),1]);
    
    for k = 1:length(z)
        for i = 1:size(X,1)
            for j = 1:size(X,2)
                Z{k}(i,j) = ls.funHan([X(i,j);Y(i,j);z(k)]);
            end
        end
    end

    % construct polygons for the both sides of the equality constraint
    pgon1 = cell(length(z),1); pgon2 = cell(length(z),1);
    contAll = cell(length(z),1);

    for k = 1:length(z)
        
        % compute contour line
        M = contourc(x,y,Z{k},[0,0]); 

        % split into connected segments
        cont = {}; cnt = 1;

        while ~isempty(M)
            len = M(2,1); cont{end+1} = [M(:,2:1+len);cnt*ones(1,len)]; 
            M = M(:,2+len:end); cnt = cnt + 1;
        end

        % compute polygons for both sides for each segment
        for i = 1:length(cont)

            % compute polygons
            [pgon1_,pgon2_] = aux_getPolygons(cont{i}(1:2,:),ls,z(k),lim);

            % combine with previous polygons
            if i == 1
                pgon1{k} = pgon1_; pgon2{k} = pgon2_; 
            else
                if isIntersecting(pgon1{k},pgon1_)
                    pgon1{k} = pgon1{k} & pgon1_;
                else
                    pgon1{k} = pgon1{k} | pgon1_;
                end
                if isIntersecting(pgon2{k},pgon2_)
                    pgon2{k} = pgon2{k} & pgon2_;
                else
                    pgon2{k} = pgon2{k} | pgon2_;
                end
            end
        end

        % store contour
        contAll{k} = [cont{:}];
    end

    % construct surface for each level
    F = repmat({{}},[length(z)-1,1]);

    for i = 1:length(pgon1)-1

        % check if polygons are empty
        if (representsa(pgon1{i+1},'emptySet') && representsa(pgon2{i+1},'emptySet')) || ...
            (representsa(pgon1{i},'emptySet') && representsa(pgon2{i},'emptySet'))
            continue;
        end

        % compute intersection of the polygons for the two levels
        p = (pgon1{i} & pgon2{i+1}) | (pgon1{i+1} & pgon2{i});

        % compute triangulation of the intersection
        if ~representsa(p,'emptySet')
            T = triangulation(p);
        else
            T = [];
        end

        % construct surface for each triangle
        for j = 1:length(T)

            % get vertices of the triangle
            V = T{j}.set.Vertices'; V = [V;zeros(1,size(V,2))];

            % add the corresponding z-coordinate to each vertex
            for k = 1:size(V,2)
                if ismember(V(1:2,k)',contAll{i+1}(1:2,:)','rows')
                    V(end,k) = z(i+1);
                else
                    V(end,k) = z(i);
                end
            end

            % store the resulting face
            F{i}{end+1} = V;
        end

        % determine contour points that are not part of the intersection
        % (because they lie directly above each other)
        ind1 = find(contains(p,contAll{i}(1:2,:),'exact',1e-5) == 0);
        ind2 = find(contains(p,contAll{i+1}(1:2,:),'exact',1e-5) == 0);

        if ~isempty(ind1) && ~isempty(ind2)

            for k = 1:length(ind1)-1
    
                % check if points are neighbours and belong to same contour
                if ind1(k+1) == ind1(k) + 1 && ...
                        contAll{i}(3,ind1(k)) == contAll{i}(3,ind1(k+1))
                    
                    % determine closed point in upper level for each point
                    [~,i1] = min(sum((contAll{i+1}(1:2,ind2) - ...
                                             contAll{i}(1:2,ind1(k))).^2));
                    [~,i2] = min(sum((contAll{i+1}(1:2,ind2) - ...
                                           contAll{i}(1:2,ind1(k+1))).^2));
    
    
                    % check if point in upper level are neighbours
                    if abs(ind2(i1) - ind2(i2)) == 1 && ...
                       contAll{i+1}(3,ind2(i1)) == contAll{i+1}(3,ind2(i2))
    
                        p1 = contAll{i}(1:2,ind1(k));
                        p2 = contAll{i}(1:2,ind1(k+1));
                        p3 = contAll{i+1}(1:2,ind2(i1));
                        p4 = contAll{i+1}(1:2,ind2(i2));
                        
                        F{i}{end+1} = [p1,p2,p3;z(i),z(i),z(i+1)];
                        F{i}{end+1} = [p3,p4,p2;z(i+1),z(i+1),z(i)];
                    end
                end
            end
        end
    end
end

function F = aux_getFaces3Dsolvable(ls,lim,ind)
% compute faces that represent an equality constraint in 3D for the case
% that the equalilty constraint is solvable for the third variable

    % setup meshgrid for x-y-plane
    N = 50; 

    ind = [setdiff(1:3,ind),ind];
    
    dx = (lim(ind(1),2)-lim(ind(1),1))/N;
    dy = (lim(ind(2),2)-lim(ind(2),1))/N;
    dz = (lim(ind(3),2)-lim(ind(3),1))/N;
    
    x = lim(ind(1),1):dx:lim(ind(1),2);
    y = lim(ind(2),1):dy:lim(ind(2),2);
    z = lim(ind(3),1):dz:lim(ind(3),2);
    
    [X,Y] = meshgrid(x,y);

    % compute z-coordinates for all grid values
    Z = zeros(size(X));
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            tmp = zeros(3,1); tmp(ind(1:2)) = [X(i,j);Y(i,j)];
            Z(i,j) = ls.solved{ind(end)}.funHan{1}.eq(tmp);
        end
    end

    % replace numbers outside of the range with NaN to prevent them from
    % getting plottet
    Z_ = inf*ones(size(Z,1)+2,size(Z,2)+2);
    Z_(2:end-1,2:end-1) = Z;

    for i = 2:size(Z_,1)-1
        for j = 2:size(Z_,2)-1
            tmp = Z_([i-1,i,i+1],[j-1,j,j+1]);
            if all(all(tmp < lim(3,1) | tmp > lim(3,2)))
                Z(i-1,j-1) = NaN;
            end
        end
    end

    % loop over all grid cells
    F_ = {};

    for i = 1:size(X,1)-1
        for j = 1:size(X,2)-1
        
            % first triangle
            if all(~isnan([Z(i,j),Z(i+1,j),Z(i,j+1)]))
                F_{end+1} = [X(i,j),X(i+1,j),X(i,j+1); ...
                             Y(i,j),Y(i+1,j),Y(i,j+1); ...
                             Z(i,j),Z(i+1,j),Z(i,j+1)];
            end

            % second triangle
            if all(~isnan([Z(i+1,j),Z(i,j+1),Z(i+1,j+1)]))
                F_{end+1} = [X(i+1,j),X(i,j+1),X(i+1,j+1); ...
                             Y(i+1,j),Y(i,j+1),Y(i+1,j+1); ...
                             Z(i+1,j),Z(i,j+1),Z(i+1,j+1)];
            end
        end
    end

    % add the faces to the correct levels
    F = repmat({{}},[length(z)-1,1]); 
    index = find(ind == 3);

    for i = 1:length(F)
        for j = 1:length(F_)
            if aux_isIntersecting1D(z(i),z(i+1),min(F_{j}(index,:)), ...
                                                  max(F_{j}(index,:)),eps)
                tmp = zeros(size(F_{j}));
                tmp(ind,:) = F_{j};
                F{i}{end+1} = tmp;
            end
        end
    end
end

function res = aux_isIntersecting1D(inf1,sup1,inf2,sup2,tol)
% check if two one-dimensional intervals intersect
    res = false;

    if inf1 <= inf2 || withinTol(inf1,inf2,tol)
        if inf2 <= sup1 || withinTol(inf2,sup1,tol)
            res = true;
        end
        
    else % inf2 < inf1
        if inf1 <= sup2 || withinTol(inf1,sup2,tol)
            res = true;
        end
    end
end

function [pgon1,pgon2] = aux_getPolygons(cont,ls,z,lim)
% compute polygons for both sides for each segment

    % construct rectangle for the domain
    V = [lim(1,1),lim(1,1),lim(1,2),lim(1,2); ...
         lim(2,1),lim(2,2),lim(2,2),lim(2,1)];

    R = polygon(V(1,:),V(2,:));

    % check if the contour is empty
    if isempty(cont)
        if ls.funHan(center(R)) <= 0
            pgon1 = R; pgon2 = []; return;
        else
            pgon1 = []; pgon2= []; return;
        end
    end

    % map start and end point of the segment to an edge of the rectangle
    e1 = aux_mapToEdge(cont(:,1),lim);
    e2 = aux_mapToEdge(cont(:,end),lim);

    % construct polygon for the first side
    if e1 == 0 || e2 == 0       % closed polygon
        pgon1 = polygon(cont(1,:),cont(2,:));
    elseif e1 == e2             % polygon that intersects one boundary
        pgon1 = polygon(cont(1,:),cont(2,:));
    else                        % polygon that intersects the boundaries
        if e2 > e1
            tmp = [cont,fliplr(V(:,e1:e2-1))];
        else
            tmp = [cont,fliplr([V(:,e1:end),V(:,1:e2-1)])];
        end
        pgon1 = polygon(tmp(1,:),tmp(2,:));
    end

    % check if first polygon belongs to the inside or outside
    c = center(pgon1); inside = false;

    if ~isempty(z)
        p = [c;z];
    else
        p = c;
    end

    if size(cont,2) == 1            % only a single point
        g = ls.der.grad(p);
        e = aux_mapToEdge(c,lim);
        if e == 0                       % point located in the interior
            inside = true;
        else                            % point located on the boundary
            d = center(R) - c;
            if g(1:2)'*d > 0
                inside = true;
            end
        end
    elseif all(min(abs(repmat(cont,[4,1]) - ...  % all points on boundary
                reshape(V,[numel(V),1])),[],1) == 0)
        g = ls.der.grad(p);
        d = center(R) - c;
        if g(1:2)'*d > 0
            inside = true;
        end
    else                            % full-dimensional polygon

        V_ = pgon1.set.Vertices'; V_ = [V_,V_(:,1)]; direc = [];

        for i = 1:size(V_,2)-1
            if ismember(V_(:,i)',cont','rows') && ...
                                    ismember(V_(:,i+1)',cont','rows')

                % check if the gradient points to the inside
                d = V_(:,i+1)-V_(:,i); d = [d(2);-d(1)];
                p1 = V_(:,i); p2 = V_(:,i+1);
                if ~isempty(z)
                    p1 = [p1;z]; p2 = [p2;z];
                end
                g1 = ls.der.grad(p1); g2 = ls.der.grad(p2);
                direc = [direc,sign(g1(1:2)'*d),sign(g2(1:2)'*d)];
            end
        end

        if sum(direc) < length(direc)/2
            inside = true;
        end
    end

    % compute second polygon
    if inside
        pgon2 = subtract(R,pgon1);
    else
        pgon2 = pgon1;
        pgon1 = subtract(R,pgon2);
    end
end

function edge = aux_mapToEdge(x,lim)
% map the point to a specific edge of the rectangle that represents the
% domain of valid values

    edge = 0;

    if abs(x(2) - lim(2,1)) < eps
        edge = 1;
    elseif abs(x(1) - lim(1,1)) < eps
        edge = 2;
    elseif abs(x(2) - lim(2,2)) < eps
        edge = 3;
    elseif abs(x(1) - lim(1,2)) < eps
        edge = 4;
    end
end

function lsFac = aux_factorLevelSet(ls)
% try to represent the level set equation as a product of terms, because 
% p(x)*q(x) = 0 is equivalent to the union of the two level sets p(x) = 0
% and q(x) = 0

    % compute terms of the product representation
    fac = factor(ls.eq);

    % create a new level set for each term of the product
    lsFac = {};

    for i = 1:length(fac)
        if ~isempty(symvar(fac(i))) && ~ismember(fac(i),fac(1:i-1))
            lsFac{end+1} = levelSet(fac(i),ls.vars,ls.compOp);
        end
    end
end

function [line,face] = aux_intersectTriangularSurfaces(V1,V2)
% compute the line that represents the intersection of two 3D surfaces
% that are spanned by 3 points

    line = []; face = [];

    if any(any(isnan(V1))) || any(any(isinf(V1))) || ...
                        any(any(isnan(V2))) || any(any(isinf(V2)))
        return;
    end

    % compute hyperplanes for the two surfaces
    [c1,d1] = aux_points2hyperplane(V1);
    [c2,d2] = aux_points2hyperplane(V2);
    
    % check if the surfaces cross the hyperplanes
    tmp1 = c2'*V1 - d2;
    
    if min(tmp1) > 0 || max(tmp1) < 0
        return;
    end
    
    tmp2 = c1'*V2 - d1;
    
    if min(tmp2) > 0 || max(tmp2) < 0
        return;
    end

    % catch the special case where the triangle lie in the same plane
    if all(abs(tmp1) < eps) && all(abs(tmp2) < eps)

        A = gramSchmidt(c1); b = V1(:,1);
        V1_ = A'*(V1 - b); poly1 = polygon(V1_(2,:),V1_(3,:));
        V2_ = A'*(V2 - b); poly2 = polygon(V2_(2,:),V2_(3,:));
        
        if ~isIntersecting(poly1,poly2)
            return;
        end

        V = vertices(poly1 & poly2);
        face = A*[zeros(1,size(V,2));V] + b;
        return;
    end
    
    % compute intersection of edges of triangle 1 with the hyperplane 2
    V1_ = [V1,V1(:,1)]; tmp1 = [tmp1,tmp1(1)]; p1 = [];
    
    for i = 1:size(V1_,2)-1
        if sign(tmp1(i)) ~= sign(tmp1(i+1))
            
            % compute intersection
            n = V1_(:,i+1) - V1_(:,i);
            
            lambda = (d2 - c2'*V1_(:,i))/(c2'*n);
            p_ = V1_(:,i) + lambda * n; 
    
            % check if point belongs to triangle 1
            if lambda >= 0 - eps && lambda <= 1 + eps
                p1 = [p1,p_];
            end
        end
    end
    
    % compute intersection of edges of triangle 2 with the hyperplane 1
    V2_ = [V2,V2(:,1)]; tmp2 = [tmp2,tmp2(1)]; p2 = [];
    
    for i = 1:size(V2_,2)-1
        if sign(tmp2(i)) ~= sign(tmp2(i+1))
            
            % compute intersection
            n = V2_(:,i+1) - V2_(:,i);
            
            lambda = (d1 - c1'*V2_(:,i))/(c1'*n);
            p_ = V2_(:,i) + lambda * n; 
    
            % check if point belongs to triangle 2
            if lambda >= 0 - eps && lambda <= 1 + eps
                p2 = [p2,p_];
            end
        end
    end
    
    % determine points with maximum distance
    p = [p1,p2]; distMax = -inf;
    
    for i = 1:size(p,2)
        for j = i+1:size(p,2)
            dist_ = sqrt(sum((p(:,i) - p(:,j)).^2));
            if dist_ > distMax
                distMax = dist_; ind = i;
            end
        end
    end
    
    dist = sqrt(sum((p - p(:,ind)).^2))/distMax;
    
    % compute intersection of ranges on the line passing through both triangles
    lb1 = min(dist(1:size(p1,2)));
    ub1 = max(dist(1:size(p1,2)));
    lb2 = min(dist(size(p1,2)+1:end));
    ub2 = max(dist(size(p1,2)+1:end));
    
    lb = max(lb1,lb2); ub = min(ub1,ub2);
    
    if ub < lb
        return;
    end
    
    [~,ind1] = min(abs(dist - lb));
    [~,ind2] = min(abs(dist - ub));
    
    line = p(:,[ind1,ind2]);
end

function [c,d] = aux_points2hyperplane(X)
% construct a hyperplane that runs through all points

    c = ndimCross([X(:,1)-X(:,2),X(:,1)-X(:,3)]);
    c = c/norm(c);
    d = c'*X(:,1);
end

function seg = aux_removeSegmentsIneqCons(cont,ineq)
% remove all points from the lines that violate the given inequality
% constraints

    seg = {};

    for i = 1:length(cont)

        inside = contains(ineq,cont{i});

        ind = find(inside(2:end) ~= inside(1:end-1));
        ind = [0,ind,length(inside)];
        
        for j = 1:length(ind)-1
            if all(inside(ind(j)+1:ind(j+1)))
                seg{end+1} = cont{i}(:,ind(j)+1:ind(j+1));
            end
        end
    end
end

function F = aux_intersect3Dinequality(F,ineq)
% intersect a triangular surface in 3D with a bunch of inequality
% constraints

    % loop over all edges of the triangular surface and determine points
    % where the inequality constraints become true
    V = [F,F(:,1)]; m = [];
    inside = contains(ineq,V);

    for i = 1:3
        if inside(i) ~= inside(i+1)
            p = V(:,i) + (0:0.01:1).*(V(:,i+1)-V(:,i));
            tmp = ineq.funHan(p);
            ind = find(max(tmp,[],1) <= 0);
            if ~isempty(ind)
                [~,index] = min(abs(tmp(:,ind)));
                m = [m,p(:,ind(index))];
            end
        end
    end

    % construct new triangular surfaces
    ind = find(inside(1:3));

    if length(ind) == 1
        F = {[V(:,ind),m]};
    else
        F = {[V(:,ind(1)),m];[V(:,ind(2)),m]};
    end
end

function ls = aux_selectCorrectDimensions(ls,dims)
% modify the level set such that the dimension that should be plotted are
% the first ones and all other dimensions are replaced by zeros

    if length(ls.vars) ~= length(dims) || any(dims ~= 1:length(ls.vars))
        x = sym('x',[length(dims),1]);
        x_ = sym(zeros(length(ls.vars),1));
        x_(dims) = x;
    
        ls = levelSet(ls.funHan(x_),x_(x_ ~= 0),ls.compOp);
    end
end

function F = aux_combineFaces(F1,F2)
% combine two lists storing faces for each level

    if isempty(F1)
        F = F2; return;
    end

    if isempty(F2)
        F = F1; return;
    end

    F = F1;

    for i = 1:length(F)
        F{i} = [F{i},F2{i}];
    end
end

% ------------------------------ END OF CODE ------------------------------
