function rotation = s2_offset(points_1)
%S2_OFFSET Experimental offset rotation of S^2
%
%Syntax
% rotation = s2_offset(points_1);
%
%Description
% ROTATION = S2_OFFSET(POINTS_1) sets ROTATION to be an R^3 rotation matrix which
% rotates the north pole of S^2 to a point specified by the points of POINTS_1.
%
% POINTS_1 must be a 2 by M matrix, representing M points of S^2 in spherical
% polar coordinates, with M a positive integer.
%
%Examples
% > s
% s =
%          0    0.7854    2.3562    3.9270    5.4978         0
%          0    1.5708    1.5708    1.5708    1.5708    3.1416
%
% > r=s2_offset(s)
% r =
%     0.0000    0.7071    0.7071
%    -1.0000    0.0000         0
%    -0.0000   -0.7071    0.7071
%
%See also
% ROT3, CIRCLE_OFFSET

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_in_collar = size(points_1,2);
if n_in_collar > 2
    if (n_in_collar > 3) && (points_1(2,2) == points_1(2,3))
        a_3 = (points_1(1,2) + points_1(1,3))/2;
    else
        a_3 = points_1(1,2) + pi;
    end
    a_2 = points_1(2,2)/2;
else
    a_3 = 0;
    a_2 = pi/2;
end
rotation = rot3(2,-a_2)*rot3(3,-a_3);
%
%end function
