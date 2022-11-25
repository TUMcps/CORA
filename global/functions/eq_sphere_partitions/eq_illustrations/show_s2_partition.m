function [movie_frame] = show_s2_partition(N,varargin)
%SHOW_S2_PARTITION 3D illustration of an EQ partition of S^2
%
%Syntax
% [movie_frame] = show_s2_partition(N,options);
%
%Description
% SHOW_S2_PARTITION(N) uses a 3d plot to illustrate the partition of
% the unit sphere S^2 into N regions.
%
% MOVIE_FRAME = SHOW_S2_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% SHOW_S2_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% SHOW_S2_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% SHOW_S2_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% SHOW_S2_PARTITION(N,'title','show')
% SHOW_S2_PARTITION(N,'title','hide')
% Show or hide title (default 'show').
%
% SHOW_S2_PARTITION(N,'points','show')
% SHOW_S2_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% SHOW_S2_PARTITION(N,'sphere','show')
% SHOW_S2_PARTITION(N,'sphere','hide')
% Show or hide the unit sphere S^2 (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
% > show_s2_partition(10)
% > frames=show_s2_partition(9,'offset','extra')
% frames =
% 1x10 struct array with fields:
%     cdata
%     colormap
% > show_s2_partition(99,'points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, PROJECT_S2_PARTITION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;
popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 16;
gdefault.show_title  = false;
gdefault.show_points = true;
gdefault.show_sphere = true;
gopt = illustration_options(gdefault, varargin{:});

dim = 2;

surf_jet;

if gopt.show_title
    if gopt.show_points
        pointstr = ', showing the center point of each region';
    else
        pointstr = '';
    end
    titlestr = sprintf(...
        '\nRecursive zonal equal area partition of {S^2} \n into %d regions%s.',...
        N,pointstr);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

frame_no = 1;
if nargout > 0
    movie_frame(frame_no) = getframe(gcf);
    frame_no = frame_no + 1;
end

if gopt.show_sphere
    show_s2_sphere;
    hold on
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
        frame_no = frame_no + 1;
    end
end

R = eq_regions(dim,N,popt.extra_offset);
top_colat = 0;
for i = N:-1:2
    if top_colat ~= R(2,1,i)
        top_colat = R(2,1,i);
        pause(0);
    end
    show_s2_region(R(:,:,i),N);
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
        frame_no = frame_no + 1;
    end
end

if gopt.show_points
    x = eq_point_set(dim,N,popt.extra_offset);
    show_r3_point_set(x,'sphere','hide','title','hide');
    hold on
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
        frame_no = frame_no + 1;
    end
end

hold off
%
% end function

function show_s2_region(region,N)
%SHOW_S2_REGION Illustrate a region of S^2
%
%Syntax
% show_s2_region(region,N);
%
%Description
% SHOW_S2_REGION(REGION,N) uses 3D surface plots to illustrate a region of S^2.
% The region is given as a 2 x 2 matrix in spherical polar coordinates

tol = eps*2^5;

dim = size(region,1);
t = region(:,1);
b = region(:,2);

if abs(b(1)) < tol
    b(1) = 2*pi;
end
pseudo = 0;
if abs(t(1)) < tol && abs(b(1)-2*pi) < tol
    pseudo = 1;
end
n = 21;
delta = 1/(n-1);
h = 0:delta:1;
t_to_b = zeros(dim,n);
b_to_t = t_to_b;
r = sqrt(1/N)/12;
for k = 1:dim
    if ~pseudo || k < 2
        L = 1:dim;
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:) = t(j(1))+(b(j(1))-t(j(1)))*h;
        t_to_b(j(2),:) = t(j(2))*ones(1,n);
        t_to_b_x = polar2cart(t_to_b);
        %plot3(t_to_b_x(1,:),t_to_b_x(2,:),t_to_b_x(3,:),'LineWidth',1,'Color','k');
        [X,Y,Z] = fatcurve(t_to_b_x,r);
        surface(X,Y,Z,-ones(size(Z)),...
       'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        axis equal
        hold on
    end
end
grid off
axis off
%
% end function
