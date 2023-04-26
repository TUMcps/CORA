function project_point_set(points,varargin)
%PROJECT_POINT_SET Use projection to illustrate a point set of S^2 or S^3
%
%Syntax
% project_point_set(points,options);
%
%Description
% PROJECT_POINT_SET(POINTS,OPTIONS) uses projection to illustrate a point set
% of S^2 or S^3, represented by the Cartesian coordinates POINTS.
%
% The argument POINTS must be an array of real numbers of size (3 by N) or
% (4 by N), where N is a positive integer.
%
% PROJECT_POINT_SET(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_POINT_SET(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% PROJECT_POINT_SET(N,'title','show')
% PROJECT_POINT_SET(N,'title','hide')
% Show or hide titles (default 'hide').
%
% PROJECT_POINT_SET(N,'proj','stereo')
% PROJECT_POINT_SET(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
%Notes
% The points are assumed to all lie on the unit sphere S^dim, where dim == 2 or
% dim == 3. The first point POINTS(:,1) should be the North pole [1,0,0]' or 
% [1,0,0,0]'.
%
%Examples
% > x
%  x =
%           0    0.0000   -0.0000    0.0000
%           0    1.0000   -1.0000         0
%      1.0000    0.0000    0.0000   -1.0000
%  
% > project_point_set(x)
%
%See also
% ILLUSTRATION_OPTIONS

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

gdefault.fontsize = 16;
gdefault.show_title  = false;
gdefault.stereo = true;

gopt = illustration_options(gdefault, varargin{:});

if gopt.stereo
    projection = 'x2stereo';
else
    projection = 'x2eqarea';
end

dim = size(points,1)-1;
N = size(points,2);

switch dim
case 2
    t = feval(projection,points);
    r = real(pi-acos(points(end,N-size(t,2)+1:end)));
    limit = pi;
    if N <= 4
        c = 'k';
    else
        c = r;
    end
    s = ceil(40*r/limit)+8;
    h = scatter(t(1,:),t(2,:),s,c,'filled');
    if N < 256
        set(h,'MarkerEdgeColor',[0.5,0.5,0.5]);
    else
        set(h,'MarkerEdgeColor','none');
    end
    axis equal
    grid off
    axis off
case 3
    t = feval(projection,points);
    r = real(pi-acos(points(end,N-size(t,2)+1:end)));
    if gopt.stereo
        limit = pi;
    else
        limit = 2*pi;
    end
    s = (r+1)/(limit*10);
    sphere_detail = ceil(20/log(N+1));
    [X,Y,Z] = sphere;

    for k = 1:size(t,2)
        surf(t(1,k)+s(k)*X,t(2,k)+s(k)*Y,t(3,k)+s(k)*Z,r(k)*ones(size(X)),...
             'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        axis equal
        hold on
    end
    grid off
    axis off
    camlight right
otherwise
    error('project_point_set(points_x): points_x must be a point set in R^3 or or R^4');
end

if gopt.show_title
    titlestr = sprintf(...
        'Point set containing %d points.',N);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

hold off
%
% end function
