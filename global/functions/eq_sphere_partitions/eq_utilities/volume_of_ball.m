function volume = volume_of_ball(dim)
%VOLUME_OF_BALL Volume of the unit ball
%
%Syntax
% volume = volume_of_ball(dim);
%
%Description
% VOLUME = VOLUME_OF_BALL(dim) sets VOLUME to be the volume of the unit ball
% B^dim in R^dim which is enclosed by the sphere S^(dim-1).
%
% The argument dim must be a positive integer or an array of positive integers.
% The result VOLUME will be an array of the same size as dim.
%
%Notes
% The volume of B^dim is defined via the Lebesgue measure on R^dim.
%
% Ref: [WeiMW].
%
%Examples
% > volume=volume_of_ball(1:7)
% volume =
%     2.0000    3.1416    4.1888    4.9348    5.2638    5.1677    4.7248
%
%See also
% AREA_OF_SPHERE

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

volume = area_of_sphere(dim-1)./dim;
%
% end function
