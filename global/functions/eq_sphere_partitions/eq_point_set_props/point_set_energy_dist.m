function [energy,min_dist] = point_set_energy_dist(points,s)
%POINT_SET_ENERGY_DIST Energy and minimum distance of a point set
%
%Syntax
% [energy,min_dist] = point_set_energy_dist(points,s);
%
%Description
% [ENERGY,MIN_DIST] = POINT_SET_ENERGY_DIST(POINTS,s) sets ENERGY to be the
% energy of the r^(-s) potential on the point set POINTS, and sets MIN_DIST 
% to be the minimum Euclidean distance between points of POINTS.
%
% POINTS must be an array of real numbers of size (M by N), where M and N
% are positive integers, with each of the N columns representing a point of
% R^M in Cartesian coordinates.
% The result MIN_DIST is optional.
%
% [ENERGY,MIN_DIST] = POINT_SET_ENERGY_DIST(POINTS) uses the default value of
% dim-1 for s.
%
%Notes
% The value of ENERGY for a single point is 0.
% Since this function is usually meant to be used for points on a unit sphere,
% the value of MIN_DIST for a single point is defined to be 2.
%
%Examples
% > x
% x =
%          0         0    0.0000         0
%          0         0         0         0
%          0    1.0000   -1.0000    0.0000
%     1.0000    0.0000    0.0000   -1.0000
%
% > [e,d]=point_set_energy_dist(x)
% e =
%     2.5000
% d =
%     1.4142
%
%See also
% EUCLIDEAN_DIST, EQ_ENERGY_DIST, POINT_SET_MIN_DIST

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-22 $
% Documentation files renamed
% Fix nasty but obvious bug by using separate variables dist and min_dist
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(points,1)-1;
%
% The default value of s is dim-1.
%
if nargin < 2
    s = dim-1;
end
energy = 0;
if nargout > 1
    min_dist = 2;
end
n_points = size(points,2);
for i = 1:(n_points-1)
    j = (i+1):n_points;
    dist = euclidean_dist(points(:,i)*ones(1,n_points-i),points(:,j));
    energy = energy + sum(potential(s,dist));
    if nargout > 1
        min_dist = min(min_dist,min(dist));
    end
end
%
% end function

function pot = potential(s,r)
%POTENTIAL r^(-s) potential at Euclidean radius r.
%
%Syntax
% pot = potential(s,r);

switch s
    case 0
        pot = -log(r);
    otherwise
        pot = r.^(-s);
end
%
% end function
