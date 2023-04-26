function min_dist = point_set_min_dist(points)
%POINT_SET_MIN_DIST Minimum distance between points of a point set 
%
%Syntax
% min_dist = point_set_min_dist(points);
%
%Description
% MIN_DIST = POINT_SET_MIN_DIST(POINTS) sets MIN_DIST to be the minimum
% Euclidean distance between points of the point set POINTS.
%
% POINTS must be an array of real numbers of size (M by N), where M and N
% are positive integers, with each of the N columns representing a point of
% R^M in Cartesian coordinates.
%
%Notes
% Since this function is usually meant to be used for points on a unit sphere, 
% the value of POINT_SET_MIN_DIST for a single point is defined to be 2.
%
%Examples
% > x
% x =
%          0         0    0.0000         0
%          0         0         0         0
%          0    1.0000   -1.0000    0.0000
%     1.0000    0.0000    0.0000   -1.0000
%  
% > d=point_set_min_dist(x)
% d =
%     1.4142
%
%See also
% EUCLIDEAN_DIST, MIN_EUCLIDEAN_DIST, POINT_SET_ENERGY_DIST

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

min_dist = 2;
n_points = size(points,2);
for i = 1:(n_points-1)
    j = (i+1):n_points;
    dist = euclidean_dist(points(:,i)*ones(1,n_points-i),points(:,j));
    min_dist = min(min_dist,min(dist));
end
%
% end function
