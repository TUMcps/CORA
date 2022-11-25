function area = area_of_ideal_region(dim,N)
%AREA_OF_IDEAL_REGION Area of one region of an EQ partition
%
%Syntax
% area = area_of_ideal_region(dim,N);
%
%Description
% AREA = AREA_OF_IDEAL_REGION(dim,N) sets AREA to be the area of one of N equal
% area regions on S^dim, that is 1/N times AREA_OF_SPHERE(dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as N.
%
%Examples
% > area=area_of_ideal_region(3,1:6)
% area =
%    19.7392    9.8696    6.5797    4.9348    3.9478    3.2899
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

area = area_of_sphere(dim)./N;
%
% end function
