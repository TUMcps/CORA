function offset = circle_offset(n_top,n_bot,extra_twist)
%CIRCLE_OFFSET Try to maximize minimum distance of center points for S^2 collars
%
% Given n_top and n_bot, calculate an offset.
%
% The values n_top and n_bot represent the numbers of 
% equally spaced points on two overlapping circles.
% The offset is given in multiples of whole rotations, and
% consists of three parts;
% 1) Half the difference between a twist of one sector on each of bottom and top. 
% This brings the centre points into alignment.
% 2) A rotation which will maximize the minimum angle between
% points on the two circles.
% 3) An optional extra twist by a whole number of sectors on the second circle.
% The extra twist is added so that the location of
% the minimum angle  between circles will
% progressively twist around the sphere with each collar.
%
% offset = circle_offset(n_top,n_bot,extra_twist);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if nargin < 3
    extra_twist = false;
end
offset = (1/n_bot - 1/n_top)/2 + gcd(n_top,n_bot)/(2*n_top*n_bot);
if extra_twist
    twist = 6;
    offset = offset + twist/n_bot;
end
% end function
