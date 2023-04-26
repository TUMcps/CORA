function result = x2eqarea(x)
%X2EQUAREA Equal area projection of 3D Euclidean points
%
% result = x2eqarea(x);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(x,1)-1;
theta = acos(-x(dim+1,:));
r = (area_of_cap(dim,theta)/volume_of_ball(dim)).^(1/dim);
mask = (sin(theta) == 0);
scale = ones(dim,1)*(r(~mask)./sin(theta(~mask)));
eqarea = zeros(dim,size(x,2));
eqarea(:,~mask) = x(1:dim,~mask).*scale;
mask1 = (x(dim+1,:) == 1);
result = eqarea(:,~mask1);
% end function
