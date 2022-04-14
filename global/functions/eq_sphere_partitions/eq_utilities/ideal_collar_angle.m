function angle = ideal_collar_angle(dim,N)
%IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
%
%Syntax
% angle = ideal_collar_angle(dim,N);
%
%Description
% ANGLE = IDEAL_COLLAR_ANGLE(dim,N) sets ANGLE to the ideal angle for the 
% spherical collars of an EQ partition of the unit sphere S^dim into N regions.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers. 
% The result ANGLE will be an array of the same size as N.
%
%Notes
% The ideal collar angle is determined by the side of a dim-dimensional
% hypercube of the same volume as the area of a single region of an N region
% equal area partition of S^dim.
%
% Since the EQ partition for N < 3 has no spherical collars, 
% the recursive zonal equal area sphere partitioning algorithm does not use 
% ideal_collar_angle(dim,N) for N < 3.
%
%Examples
% > angle = ideal_collar_angle(2,10)
%  angle =
%      1.1210
%  
% > angle = ideal_collar_angle(3,1:6)
%  angle =
%      2.7026    2.1450    1.8739    1.7025    1.5805    1.4873
%
%See also
% AREA_OF_IDEAL_REGION

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

angle = area_of_ideal_region(dim,N).^(1/dim);
%
% end function
