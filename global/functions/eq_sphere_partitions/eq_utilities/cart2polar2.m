function s = cart2polar2(x)
%CART2POLAR2 Convert from Cartesian to spherical coordinates on sphere S^2
%
%Syntax
% s = cart2polar2(x);
%
%Description
% S = CART2POLAR2(X) sets S to be the spherical polar coordinates of the points
% represented by the Cartesian coordinates X:
% S = [phi;theta]: phi in [0, 2*pi), theta in [0, pi].
%
% The argument X must be an array of real numbers of size (3 by N), where N is 
% any positive integer. The result S will be an array of size (2 by N).
%
%Examples
% > x
% x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% > s=cart2polar2(x)
% s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
%Note
% CART2POLAR2(X) projects any X in R^3 onto the sphere S^2 via a line through 
% the origin. The origin [0 0 0]' is itself projected onto a point on the
% equator such that
%
%     POLAR2CART(CART2POLAR2([0 0 0]')) == [1 0 0]'.
%
%See also
% POLAR2CART

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.10 $ $Date 2005-05-28 $
% Use cart2sph
% Change name from x2s2 to cart2polar2
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.
 
[phi, theta] = cart2sph(x(1,:),x(2,:),x(3,:));
s = [mod(phi, 2*pi); pi/2-theta];
