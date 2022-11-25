function [s_cap,n_regions] = eq_caps(dim,N)
%EQ_CAPS Partition a sphere into to nested spherical caps
%
%Syntax
% [s_cap,n_regions] = eq_caps(dim,N);
%
%Description
% [S_CAP,N_REGIONS] = EQ_CAPS(dim,N) does the following:
% 1) partitions the unit sphere S^dim into a list of spherical caps of
%    increasing colatitude and thus increasing area,
% 2) sets S_CAP to be an array of size (1 by N_COLLARS+2),
%    containing increasing colatitudes of caps, and
% 3) sets N_REGIONS to be an array of size (1 by N_COLLARS+2),
%    containing the intger number of regions in each corresponding zone of 
%    S^dim.
%
% The argument N is assumed to be a positive integer.
%
%Notes
% The value N_COLLARS is a positive integer function of dim and N.
%
% S_CAP[1] is C_POLAR, the colatitude of the North polar cap.
% S_CAP[N_COLLARS+1] is pi-C_POLAR.
% S_CAP[N_COLLARS+2] is pi.
%
% N_REGIONS[1] is 1.
% N_REGIONS[N_COLLARS+2] is 1.
% The sum of N_REGIONS is N.
%
%Examples
% > [s_cap,n_regions] = eq_caps(2,10)
% s_cap =
%     0.6435    1.5708    2.4981    3.1416
% n_regions =
%      1     4     4     1
%  
% > [s_cap,n_regions] = eq_caps(3,6)
% s_cap =
%     0.9845    2.1571    3.1416
% n_regions =
%      1     4     1
%
%See also
% EQ_REGIONS, EQ_POINT_SET_POLAR

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
error(nargchk(2,2,nargin));
error(nargoutchk(2,2,nargout));
%
% dim is the number of dimensions
% N is the number of regions
%
if dim == 1
    %
    % We have a circle. Return the angles of N equal sectors.
    %
    sector = 1:N;
    %
    % Make dim==1 consistent with dim>1 by
    % returning the longitude of a sector enclosing the
    % cumulative sum of arc lengths given by summing n_regions.
    %
    s_cap = sector*2*pi/N;
    n_regions = ones(size(sector));
    %
elseif N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    s_cap = [pi];
    n_regions = [1];
    %
else
    %
    % Given dim and N, determine c_polar,
    % the colatitude of the North polar spherical cap.
    %
    c_polar = polar_colat(dim,N);
    %
    % Given dim and N, determine the ideal angle for spherical collars.
    % Based on N, this ideal angle, and c_polar,
    % determine n_collars, the number of collars between the polar caps.
    %
    n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
    %
    % Given dim, N, c_polar and n_collars, determine r_regions,
    % a list of the ideal real number of regions in each collar,
    % plus the polar caps.
    % The number of elements is n_collars+2.
    % r_regions[1] is 1.
    % r_regions[n_collars+2] is 1.
    % The sum of r_regions is N.
    %
    r_regions = ideal_region_list(dim,N,c_polar,n_collars);
    %
    % Given N and r_regions, determine n_regions,
    % a list of the natural number of regions in each collar and
    % the polar caps.
    % This list is as close as possible to r_regions.
    % The number of elements is n_collars+2.
    % n_regions[1] is 1.
    % n_regions[n_collars+2] is 1.
    % The sum of n_regions is N.
    %
    n_regions = round_to_naturals(N,r_regions);
    %
    % Given dim, N, c_polar and n_regions, determine s_cap,
    % an increasing list of colatitudes of spherical caps which enclose the same area
    % as that given by the cumulative sum of regions.
    % The number of elements is n_collars+2.
    % s_cap[1] is c_polar.
    % s_cap[n_collars+1] is Pi-c_polar.
    % s_cap[n_collars+2] is Pi.
    %
    s_cap = cap_colats(dim,N,c_polar,n_regions);
    %
end
%
% end function
