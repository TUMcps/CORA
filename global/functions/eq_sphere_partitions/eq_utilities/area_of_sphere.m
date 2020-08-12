function area = area_of_sphere(dim)
%AREA_OF_SPHERE Area of sphere
%
%Syntax
% area = area_of_sphere(dim);
%
%Description
% AREA = AREA_OF_SPHERE(dim) sets AREA to be the area of the sphere S^dim,
%
% The argument dim must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as dim.
%
%Notes
% The area of S^dim is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [Mue98] p39.
%
%Examples
% > area=area_of_sphere(1:7)
% area =
%     6.2832   12.5664   19.7392   26.3189   31.0063   33.0734   32.4697
%
%See also
% AREA_OF_CAP, VOLUME_OF_BALL

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

power = (dim+1)/2;
area = (2*pi.^power./gamma(power));
%
% end function
