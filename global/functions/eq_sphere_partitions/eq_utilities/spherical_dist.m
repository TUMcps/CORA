function a_dist = spherical_dist(x,y)
%SPHERICAL_DIST Spherical distance between two points on the sphere
%
%Syntax
% a_dist = spherical_dist(x,y);
%
%Description
% A_DIST = SPHERICAL_DIST(X,Y) sets A_DIST to be the spherical distance
% between the two points X and Y.
%
% The arguments X and Y must be arrays of the same size, M by N, where M and N
% are positive integers. Each of X and Y is assumed to represent N points in
% R^M, in Cartesian coordinates.
% The result A_DIST will be a 1 by N array.
%
%Examples
% > x
% x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% > y
% y =
%          0    0.0000   -0.0000    0.0000
%    -0.5000    0.8660   -0.8660    0.5000
%     0.8660    0.5000   -0.5000   -0.8660
%
% > a_dist=spherical_dist(x,y)
% a_dist =
%     0.5236    0.5236    0.5236    0.5236
%
%See also
% EUCLIDEAN_DIST, E2S, S2E

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

a_dist = acos(sum(x.*y));
%
% end function
