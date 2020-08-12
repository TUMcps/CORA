function show_r3_point_set(points_x,varargin)
%SHOW_R3_POINT_SET 3D illustration of a point set
%
%Syntax
% show_r3_point_set(POINTS_X,options);
%
%Description
% SHOW_R3_POINT_SET(POINTS_X) uses a 3d plot to illustrate a point set in relation to
% the unit sphere S^2.
%
% The argument POINTS_X must be an array of real numbers of size (3 by N), where N is a
% positive integer, representing N points of R^3.
%
% SHOW_R3_POINT_SET(POINTS_X,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% SHOW_R3_POINT_SET(POINTS_X,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% SHOW_R3_POINT_SET(POINTS_X,'title','show')
% SHOW_R3_POINT_SET(POINTS_X,'title','hide')
% Show or hide title (default 'hide').
%
% SHOW_R3_POINT_SET(POINTS_X,'sphere','show')
% SHOW_R3_POINT_SET(POINTS_X,'sphere','hide')
% Show or hide the unit sphere S^2 (default 'hide').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Note
% This function is primarily for use with the point set POINTS_X as a subset of the 
% unit sphere S^2, but this is not assumed and not checked.
% If you show the unit sphere S^2 and POINTS_X contains points closer than radius 1
% from the origin, the sphere will hide these points.
%
%Examples
% > points_x
% points_x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% > show_r3_point_set(points_x,'sphere','hide')
% > show_r3_point_set(points_x,'sphere','show')
%
%See also
% ILLUSTRATION_OPTIONS, SHOW_S2_PARTITION, PROJECT_POINT_SET

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
gdefault.show_sphere = false;

gopt = illustration_options(gdefault, varargin{:});

N = size(points_x,2);

surf_jet;

if gopt.show_title
    titlestr = sprintf(...
        '\nPoint set containing %d points.',N);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

if gopt.show_sphere
    show_s2_sphere;
    hold on
end

[X,Y,Z] = sphere;

r = min(0.05,N^(-1/2)/2);
rX = r*X; rY = r*Y; rZ = r*Z;
%plot3(points_x(1,:),points_x(2,:),points_x(3,:),'ro','MarkerSize',5);
for n = 1:N
   surf(points_x(1,n)+rX,points_x(2,n)+rY,points_x(3,n)+rZ,ones(size(rZ)),...
   'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
end
%
axis equal
axis off
grid off
hold off
