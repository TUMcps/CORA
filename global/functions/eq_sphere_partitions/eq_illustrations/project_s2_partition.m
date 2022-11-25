function [movie_frame] = project_s2_partition(N,varargin)
%PROJECT_S2_PARTITION Use projection to illustrate an EQ partition of S^2
%
%Syntax
% [movie_frame] = project_s2_partition(N,options);
%
%Description
% PROJECT_S2_PARTITION(N) uses projection to illustrate the partition of
% the unit sphere S^2 into N regions.
%
% MOVIE_FRAME = PROJECT_S2_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% PROJECT_S2_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% PROJECT_S2_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_S2_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% PROJECT_S2_PARTITION(N,'title','long')
% PROJECT_S2_PARTITION(N,'title','short')
% Use long or short titles (default 'long').
%
% PROJECT_S2_PARTITION(N,'proj','stereo')
% PROJECT_S2_PARTITION(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% PROJECT_S2_PARTITION(N,'points','show')
% PROJECT_S2_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
% > project_s2_partition(10)
% > frames=project_s2_partition(9,'offset','extra','proj','eqarea')
% frames =
% 1x10 struct array with fields:
%     cdata
%     colormap
% > project_s2_partition(99,'proj','eqarea','points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, SHOW_S2_PARTITION,
% PROJECT_S3_PARTITION

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
gdefault.stereo =        true;
gdefault.show_title =    true;
gdefault.long_title =    true;
gdefault.show_points =   true;

gopt = illustration_options(gdefault, varargin{:});

dim = 2;

Phi = 2*pi*[0:1/40:1];
X = cos(Phi);
Y = sin(Phi);
if gopt.stereo
    r = 0;
else
    r = (area_of_sphere(dim)/volume_of_ball(dim)).^(1/dim);
end
plot(r*X,r*Y,'k')
axis equal;hold on
colormap jet
grid off
axis off

if gopt.show_title
    if gopt.long_title
        if gopt.stereo
            s = 'Stereographic';
        else
            s = 'Equal area';
        end
        if gopt.show_points
            pointstr = ', showing the center point of each region';
        else
            pointstr = '';
        end

        titlestr = sprintf(...
        '%s projection of recursive zonal equal area partition of {S^2}\ninto %d regions%s.',...
        s,N,pointstr);
    else
        titlestr = sprintf('EQ(2,%d)',N);
    end
    title(titlestr, ...
    'FontWeight','bold','FontUnits','normalized','FontSize',gopt.fontsize/512);
end

if nargout > 0
    movie_frame(1) = getframe(gcf);
end

R = eq_regions(dim,N,popt.extra_offset);
for i = N:-1:2
    project_s2_region(R(:,:,i),gopt.stereo);
    if nargout > 0
        movie_frame(N-i+2) = getframe(gcf);
    end
end

if gopt.show_points
    project_s2_eq_point_set(N,popt.extra_offset,gopt.stereo);
    if nargout > 0
        movie_frame(N+1) = getframe(gcf);
    end
end

hold off
%
% end function

function project_s2_region(region, stereo)
%PROJECT_S2_REGION Use projection to illustrate an EQ region of S^2
%
%Syntax
% project_s2_region(region, stereo);
%
%Notes
% The region is given as a 2 x 2 matrix in spherical polar coordinates
%
% The default is to use stereographic projection
% If the optional second argument, stereo is false,
% then use a equal area projection.

if nargin < 2
    stereo = true;
end
if stereo
    projection = 'x2stereo';
else
    projection = 'x2eqarea';
end

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
for k = 1:dim
    if ~pseudo || k < 2
        L = 1:dim;
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:) = t(j(1))+(b(j(1))-t(j(1)))*h;
        t_to_b(j(2),:) = t(j(2))*ones(1,n);
        t_to_b_x = polar2cart(t_to_b);
        s = feval(projection,t_to_b_x);
        plot(s(1,:),s(2,:),'k');
        axis equal; hold on
        if pseudo
            axis equal; hold on
            b_to_t(j(1),:,:) = b(j(1))-(b(j(1))-t(j(1)))*h;
            b_to_t(j(2),:,:) = b(j(2))*ones(1,n);
            b_to_t_x = polar2cart(b_to_t);
            s = feval(projection,b_to_t_x);
            plot(s(1,:),s(2,:),'k');
        end
    end
end
grid off
axis off
%
% end function

function project_s2_eq_point_set(N,extra_offset,stereo)
%PROJECT_S2_EQ_POINT_SET Use projection to illustrate an EQ point set of S^2
%
%Syntax
% project_s2_eq_point_set(N,extra_offset,stereo);

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

x = eq_point_set(2,N,extra_offset);
project_point_set(x,'title','hide','proj',projection);
hold on
%
% end function
