function area = area_of_collar(dim, a_top, a_bot)
%AREA_OF_COLLAR Area of spherical collar
%
%Syntax
% area = area_of_collar(dim, a_top, a_bot);
%
%Description
% AREA = AREA_OF_COLLAR(dim, A_TOP, A_BOT) sets AREA to be the area of 
% an S^dim spherical collar specified by A_TOP, A_BOT, where
% A_TOP is top (smaller) spherical radius,
% A_BOT is bottom (larger) spherical radius.
%
% The argument dim must be a positive integer.
% The arguments A_TOP and A_BOT must be real numbers or arrays of real numbers,
% with the same array size.
% The result AREA will be an array of the same size as A_TOP.
%
%Notes
% A_TOP and A_BOT are assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
% > area=area_of_collar(2,0:2,1:3)
% area =
%     2.8884    6.0095    3.6056
%
%See also
% AREA_OF_CAP, AREA_OF_SPHERE

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

area = area_of_cap(dim, a_bot) - area_of_cap(dim, a_top);
%
% end function
