function distance = euclidean_dist(x,y)
%EUCLIDEAN_DIST Euclidean distance between two points in Cartesian coordinates
%
%Syntax
% distance = euclidean_dist(x,y);
%
%Description
% DISTANCE = EUCLIDEAN_DIST(X,Y) sets DISTANCE to be the Euclidean distance
% between the two points X and Y.
%
% The arguments X and Y must be arrays of the same size, M by N, where M and N
% are positive integers. Each of X and Y is assumed to represent N points in
% R^M, in Cartesian coordinates.
% The result DISTANCE will be a 1 by N array.
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
% > distance=euclidean_dist(x,y)
% distance =
%     0.5176    0.5176    0.5176    0.5176
%
%See also
% SPHERICAL_DIST, E2S, S2E

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

distance = sqrt(sum((x-y).^2));
%
% end function
