function illustrate_eq_algorithm(dim,N,varargin)
%ILLUSTRATE_EQ_ALGORITHM Illustrate the EQ partition algorithm
%
%Syntax
% illustrate_eq_algorithm(dim,N,options);
%
%Description
% ILLUSTRATE_EQ_ALGORITHM(dim,N) illustrates the recursive zonal equal area
% sphere partitioning algorithm, which partitions S^dim (the unit sphere in 
% dim+1 dimensional space) into N regions of equal area and small diameter.
%
% The illustration consists of four subplots:
% 1. Steps 1 and 2
% 2. Steps 3 to 5
% 3. Steps 6 and 7
% 4. Lower dimensional partitions (if dim == 2 or dim == 3)
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3. If dim > 3, extra offsets are not used.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,options) also recognizes a number of
% illustration options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'title','long')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'title','short')
% Use long or short titles (default 'short').
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'proj','stereo')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'points','show')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'points','hide')
% Show or hide center points (default 'show').
%
% See examples below.
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Notes
% The step numbers refer to the following steps of the the recursive zonal
% equal area sphere partitioning algorithm, which partition the sphere into 
% zones.
%
% 1. Determine the colatitudes of the North and South polar caps.
% 2. Determine an ideal collar angle.
% 3. Use the angle between the North and South polar caps and the ideal collar
%    angle to determine an ideal number of collars.
% 4. Use a rounding procedure to determine the actual number of collars,
%    given the ideal number of collars.
% 5. Create a list containing the ideal number of regions in each collar.
% 6. Use a rounding procedure to create a list containing the actual number of
%    regions in each collar, given the list containing the ideal number of
%    regions.
% 7. Create a list containing the colatitude of the top of each zone,
%    given the list containing the actual number of regions in each collar,
%    and the colatitudes of the polar caps.
%
%Examples
% > illustrate_eq_algorithm(3,99)
% > illustrate_eq_algorithm(3,99,'offset','extra','proj','eqarea')
% > illustrate_eq_algorithm(3,99,'proj','eqarea','points','hide')
%
%See also
% PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, SUBPLOT

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;

popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 16;
gdefault.show_title =    true;
gdefault.long_title =    false;
gdefault.stereo =        false;
gdefault.show_points =   true;

gopt = illustration_options(gdefault, varargin{:});
opt_args = option_arguments(popt,gopt);

subplot(2,2,1);axis off
illustrate_steps_1_2(dim,N,opt_args);

subplot(2,2,2);axis off
illustrate_steps_3_5(dim,N,opt_args);

subplot(2,2,3);axis off
illustrate_steps_6_7(dim,N,opt_args);

subplot(2,2,4);axis off
cla

gopt.fontsize = 32;
switch dim
case 2
    opt_args = option_arguments(popt,gopt);
    project_s2_partition(N,opt_args{:});
case 3
    opt_args = option_arguments(popt,gopt);
    [s,m] = eq_caps(dim,N);
    max_collar = min(4,size(m,2)-2);
    for k = 1:max_collar
        subn = 9+2*k-mod(k-1,2);
        subplot(4,4,subn);axis off
        project_s2_partition(m(1+k),opt_args{:});
    end
end
%
% end function

function illustrate_steps_1_2(dim,N,varargin)
% Illustrate steps 1 and 2 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_1_2(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});
h = [0:1/90:1];
% Plot a circle to represent dth coordinate of S^d
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);

k = [-1:1/20:1];
j = ones(size(k));

% Plot the bounding parallels of the polar caps
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2)
plot(sin(c_polar)*k,-cos(c_polar)*j,'r','LineWidth',2)

% Plot the North-South axis
plot(zeros(size(j)),k,'b','LineWidth',1)
% Plot the polar angle
plot(sin(c_polar)*h,cos(c_polar)*h,'b','LineWidth',2)

text(0.05,2/3,'\theta_c','Fontsize',gopt.fontsize);

% Plot the ideal collar angle
Delta_I = ideal_collar_angle(dim,N);
theta = c_polar + Delta_I;
plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2)

mid = c_polar + Delta_I/2;
text(sin(mid)*2/3,cos(mid)*2/3,'\Delta_I','Fontsize',gopt.fontsize);

% Plot an arc to indicate angles
theta = h*(c_polar + Delta_I);
plot(sin(theta)/5,cos(theta)/5,'b','LineWidth',1)

text(-0.9,-0.1,sprintf('V(\\theta_c) = V_R \n    = \\sigma(S^{%d})/%d',dim,N),...
    'Fontsize',gopt.fontsize);

caption_angle = min(mid + 2*Delta_I,pi-c_polar);
text(sin(caption_angle)/3,cos(caption_angle)/3,sprintf('\\Delta_I = V_R^{1/%d}',dim),...
    'Fontsize',gopt.fontsize);

if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 1 to 2\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end

hold off
%
% end function

function illustrate_steps_3_5(dim,N,varargin)
% Illustrate steps 3 to 5 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_3_5(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});

h = [0:1/90:1];
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
s_cap = cap_colats(dim,N,c_polar,r_regions);

k = [-1:1/20:1];
j = ones(size(k));
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2);

plot(zeros(size(j)),k,'b','LineWidth',1)

for collar_n = 0:n_collars
    zone_n = 1+collar_n;
    theta = s_cap(zone_n);
    plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2);
    theta_str = sprintf('\\theta_{F,%d}',zone_n);
    text(sin(theta)*1.1,cos(theta)*1.1,theta_str,'Fontsize',gopt.fontsize);
    if collar_n ~= 0
        plot(sin(theta)*k, cos(theta)*j,'r','LineWidth',2);
        theta_p = s_cap(collar_n);
        arc = theta_p + (theta-theta_p)*h;
        plot(sin(arc)/5,cos(arc)/5,'b','LineWidth',1);
        mid = (theta_p + theta)/2;
        text(sin(mid)/2,cos(mid)/2,'\Delta_F','Fontsize',gopt.fontsize);
        y_str = sprintf('y_{%d} = %3.1f...',collar_n,r_regions(zone_n));
        text(-sin(mid)+1/20,cos(mid)+(mid-pi)/30,y_str,'Fontsize',gopt.fontsize);
    end
end
if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 3 to 5\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end
hold off
%
% end function

function illustrate_steps_6_7(dim,N,varargin)
% Illustrate steps 6 to 7 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_6_7(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});

h = [0:1/90:1];
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
n_regions = round_to_naturals(N,r_regions);
s_cap = cap_colats(dim,N,c_polar,n_regions);

k = [-1:1/20:1];
j = ones(size(k));
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2);

plot(zeros(size(j)),k,'b','LineWidth',1)

for collar_n = 0:n_collars
    zone_n = 1+collar_n;
    theta = s_cap(zone_n);
    plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2);
    theta_str = sprintf('\\theta_{%d}',zone_n);
    text(sin(theta)*1.1,cos(theta)*1.1,theta_str,'Fontsize',gopt.fontsize);
    if collar_n ~= 0
        plot(sin(theta)*k, cos(theta)*j,'r','LineWidth',2);
        theta_p = s_cap(collar_n);
        arc = theta_p + (theta-theta_p)*h;
        plot(sin(arc)/5,cos(arc)/5,'b','LineWidth',1);
        mid = (theta_p + theta)/2;
        Delta_str = sprintf('\\Delta_{%i}',collar_n);
        text(sin(mid)/2,cos(mid)/2,Delta_str,'Fontsize',gopt.fontsize);
        m_str = sprintf('m_{%d} =%3.0f',collar_n,n_regions(zone_n));
        text(-sin(mid)+1/20,cos(mid)+(mid-pi)/30,m_str,'Fontsize',gopt.fontsize);
    end
end
if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 6 to 7\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end
hold off
%
% end function

function arg = option_arguments(popt,gopt)

k = 1;
if isfield(popt,'extra_offset')
    arg{k} = 'offset';
    if popt.extra_offset
        arg{k+1} = 'extra';
    else
        arg{k+1} = 'normal';
    end
    k = k+2;
end

if isfield(gopt,'fontsize')
    arg{k} = 'fontsize';
    arg{k+1} = gopt.fontsize;
    k = k+2;
end

if isfield(gopt,'stereo')
    arg{k} = 'proj';
    if gopt.stereo
        arg{k+1} = 'stereo';
    else
        arg{k+1} = 'eqarea';
    end
    k = k+2;
end    

if isfield(gopt,'show_title')
    arg{k} = 'title';
    if gopt.show_title
        if isfield(gopt,'long_title')
            if gopt.long_title
                arg{k+1} = 'long';
            else
                arg{k+1} = 'short';
            end
        else
            arg{k+1} = 'show';
        end
    else
        arg{k+1} = 'none';
    end
    k = k+2;
elseif isfield(gopt,'long_title')
    arg{k} = 'title';
    if gopt.long_title
        arg{k+1} = 'long';
    else
        arg{k+1} = 'short';
    end
    k = k+2;
end


if isfield(gopt,'show_points')
    arg{k} = 'points';
    if gopt.show_points
        arg{k+1} = 'show';
    else
        arg{k+1} = 'hide';
    end
    k = k+2;
end

if isfield(gopt,'show_surfaces')
    arg{k} = 'surf';
    if gopt.show_surfaces
        arg{k+1} = 'show';
    else
        arg{k+1} = 'hide';
    end
    k = k+2;
end    

