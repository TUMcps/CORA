function [movie_frame] = project_s3_partition(N,varargin)
%PROJECT_S3_PARTITION Use projection to illustrate an EQ partition of S^3
%
%Syntax
% [movie_frame] = project_s3_partition(N,options);
%
%Description
% PROJECT_S3_PARTITION(N) uses projection to illustrate the partition of
% the unit sphere S^3 into N regions.
%
% MOVIE_FRAME = PROJECT_S3_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% PROJECT_S3_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% PROJECT_S3_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_S3_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 18).
%
% PROJECT_S3_PARTITION(N,'title','long')
% PROJECT_S3_PARTITION(N,'title','short')
% Use long or short titles (default 'long').
%
% PROJECT_S3_PARTITION(N,'proj','stereo')
% PROJECT_S3_PARTITION(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% PROJECT_S3_PARTITION(N,'points','show')
% PROJECT_S3_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% PROJECT_S3_PARTITION(N,'surf','show')
% PROJECT_S3_PARTITION(N,'surf','hide')
% Show or hide surfaces of regions (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
% > project_s3_partition(10)
% > frames=project_s3_partition(9,'offset','extra','proj','eqarea')
% frames =
% 1x18 struct array with fields:
%     cdata
%     colormap
% > project_s3_partition(99,'proj','eqarea','points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, PROJECT_S2_PARTITION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;

popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 18;
gdefault.stereo =        true;
gdefault.long_title =    true;
gdefault.show_points =   true;
gdefault.show_surfaces = true;

gopt = illustration_options(gdefault, varargin{:});

dim = 3;

[X,Y,Z] = sphere(90);
if gopt.stereo
    r = 0;
else
    r = (area_of_sphere(dim)/volume_of_ball(dim)).^(1/dim);
end

hold off

if gopt.show_surfaces
    surf(r*X,r*Y,r*Z,zeros(size(Z)),...
        'FaceAlpha',1/20,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
else
    plot3(0,0,0,'w.')
end
axis equal;hold on
camlight right
colormap jet
grid off
axis off

if gopt.long_title
    if gopt.stereo
        s = 'Stereographic';
    else
        s = 'Equal volume';
    end
    if gopt.show_points
        pointstr = ', showing the center point of each region';
    else
        pointstr = '';
    end

    title(sprintf(...
        '\n%s projection of recursive zonal equal area partition of {S^3} \n into %d regions%s.',...
        s,N,pointstr),'FontSize',gopt.fontsize);
else
    title(sprintf('\nEQ(3,%d)',N),'FontSize',gopt.fontsize);
end

axis equal
grid off
axis off

pause(0);
if nargout > 0
    movie_frame(1) = getframe(gcf);
end

if gopt.stereo && (N == 1)
    return;
end

if popt.extra_offset
    [R,dim_1_rot] = eq_regions(dim,N,popt.extra_offset);
else
    R = eq_regions(dim,N);
end

for i = N:-1:2
    if popt.extra_offset
        project_s3_region(R(:,:,i),N,gopt.stereo,gopt.show_surfaces,dim_1_rot{i});
    else
        project_s3_region(R(:,:,i),N,gopt.stereo,gopt.show_surfaces);
    end
    pause(0);
    if nargout > 0
        movie_frame(N-i+2) = getframe(gcf);
    end
end

if gopt.show_points
    project_s3_eq_point_set(N,popt.extra_offset,gopt.stereo);
    if nargout > 0
        for k=1:min(N,40)
            movie_frame(N+k) = getframe(gcf);
        end
    end
end

hold off
%
% end function

function project_s3_region(region, N, stereo, show_surfaces, rot_matrix)
%PROJECT_S3_REGION Use projection to illustrate an EQ region of S^3
%Syntax
% project_s3_region(region, stereo, show_surfaces, rot_matrix);
%
%Notes
% The region is given as a 3 x 2 matrix in spherical polar coordinates
%
% The default is to use stereographic projection
% If the optional second argument, stereo is false,
% then use a equal area projection.

if nargin < 3
    stereo = true;
end
if stereo
    projection = 'x2stereo';
else
    projection = 'x2eqarea';
end
if nargin < 4
    show_surfaces = true;
end

offset_regions = (nargin >= 5);

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
n = 33;
delta = 1/(n-1);
h = 0:delta:1;
[h1, h2]  = meshgrid(h,h);
t_to_b = zeros(dim,n,n);
b_to_t = t_to_b;
r = N^(-1/3)/32;
for k = 1:dim
    if ~pseudo || k < 3
        L = 1:dim;
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:,:) = t(j(1))+(b(j(1))-t(j(1)))*h1;
        t_to_b(j(2),:,:) = t(j(2))+(b(j(2))-t(j(2)))*h2;
        t_to_b(j(3),:,:) = t(j(3))*ones(n,n);
        t_to_b_v = reshape(t_to_b,dim,n*n);
        if offset_regions
            t_to_b_x = polar2cart([cart2polar2(rot_matrix*polar2cart(t_to_b_v(1:dim-1,:)));t_to_b_v(dim,:)]);
        else
            t_to_b_x = polar2cart(t_to_b_v);
        end
        s = reshape(feval(projection,t_to_b_x),dim,n,n);
        degenerate = (norm(s(:,1,1)-s(:,1,2)) < tol);
        if ~degenerate && (~pseudo || k > 1)
            [X,Y,Z] = fatcurve(squeeze(s(:,1,:)),r);
            surface(X,Y,Z,zeros(size(Z)),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
            axis equal; hold on
        end
        if show_surfaces
            surf(squeeze(s(1,:,:)),squeeze(s(2,:,:)),squeeze(s(3,:,:)),t(3)*ones(n,n),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        end
        axis equal; hold on
        camlight right
        b_to_t(j(1),:,:) = b(j(1))-(b(j(1))-t(j(1)))*h1;
        b_to_t(j(2),:,:) = b(j(2))-(b(j(2))-t(j(2)))*h2;
        b_to_t(j(3),:,:) = b(j(3))*ones(n,n);
        b_to_t_v = reshape(b_to_t,dim,n*n);
        if offset_regions
            b_to_t_x = polar2cart([cart2polar2(rot_matrix*polar2cart(b_to_t_v(1:dim-1,:)));b_to_t_v(dim,:)]);
        else
            b_to_t_x = polar2cart(b_to_t_v);
        end
        s = reshape(feval(projection,b_to_t_x),dim,n,n);
        degenerate = (norm(s(:,1,1)-s(:,1,2)) < tol);
        if ~degenerate && (~pseudo || (k > 1 && abs(b(2)-pi) > tol))
            [X,Y,Z] = fatcurve(squeeze(s(:,1,:)),r);
            surface(X,Y,Z,zeros(size(Z)),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        end
        if show_surfaces && k < 2
            surf(squeeze(s(1,:,:)),squeeze(s(2,:,:)),squeeze(s(3,:,:)),t(3)*ones(n,n),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
            camlight right
        end
    end
end
colormap jet
grid off
axis off
%
% end function

function project_s3_eq_point_set(N,extra_offset,stereo)
%PROJECT_S3_EQ_POINT_SET Use projection to illustrate an EQ point set of S^3
%
%Syntax
% project_s3_eq_point_set(N,min_energy,stereo);

if nargin < 2
    extra_offset = true;
end
if nargin < 3
    stereo = true;
end
if stereo
    projection = 'stereo';
else
    projection = 'eqarea';
end

x = eq_point_set(3,N,extra_offset);
project_point_set(x,'title','hide','proj',projection);
%
% end function
