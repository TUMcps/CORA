function R = rot3(axis, angle)
%ROT3 R^3 rotation about a coordinate axis
%
%Syntax
% R = rot3(axis, angle);
%
%Description
% R = ROT3(AXIS, ANGLE) sets R to be the 3 by 3 rotation matrix corresponding to
% AXIS and ANGLE. Use this to create rotation matrices from Euler angles.
%
% AXIS must be 1, 2, or 3.
% ANGLE must be a real number.
%
%Examples
% > r=rot3(1,pi/6)
% r =
%     1.0000         0         0
%          0    0.8660   -0.5000
%          0    0.5000    0.8660
%
% > r=rot3(2,pi/6)
% r =
%     0.8660         0   -0.5000
%          0    1.0000         0
%     0.5000         0    0.8660
%
% > r=rot3(3,pi/6)
% r =
%     0.8660   -0.5000         0
%     0.5000    0.8660         0
%          0         0    1.0000

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

c = cos(angle);
s = sin(angle);
switch axis
case 1
    R = [...
            1, 0, 0; ...
            0, c, -s; ...
            0, s, c ];
case 2
    R = [...
            c, 0, -s; ...
            0, 1, 0; ...
            s, 0, c];
case 3
    R = [...
            c, -s, 0; ...
            s, c,  0; ...
            0, 0,  1];
end
%
%end function

