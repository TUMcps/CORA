function x = polar2cart(s)
%POLAR2CART Convert spherical polar to Cartesian coordinates
%
%Syntax
% x = polar2cart(s);
%
%Description
% X = POLAR2CART(S) sets X to be the Cartesian coordinates of the points
% represented by the spherical polar coordinates S.
%
% S is assumed to be an array of size (dim by N) representing N points of
% S^dim in spherical polar coordinates, where dim and N are positive integers.
% N will be an array of size (dim+1 by N).
%
%Examples
% > s
% s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% > x=polar2cart(s)
% x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
%See also
% CART2POLAR2

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(s,1);
n=size(s,2);
x = zeros(dim+1,n);
sinprod  = 1;
for k = dim:-1:2
    x(k+1,:) = sinprod.*cos(s(k,:));
    sinprod  = sinprod.*sin(s(k,:));
end
x(2,:)=sinprod.*sin(s(1,:));
x(1,:)=sinprod.*cos(s(1,:));
r = sqrt(sum(x.^2));
mask = (r ~= 1);
if  size(r(mask),2) > 0
    x(:,mask) = x(:,mask)./(ones(dim+1,1)*r(mask));
end
%
% end function
