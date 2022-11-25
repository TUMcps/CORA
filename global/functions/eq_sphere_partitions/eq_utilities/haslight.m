function lit = haslight(axish)
%HASLIGHT Check if axis handle has a light attached
%
%Syntax
% lit = haslight(axish);
%
%Description
% LIT = HASLIGHT(AXISH) sets LIT to true if the axis specified by axis
% handle AXISH has a child of type 'light', and sets LIT to false otherwise.
%
%Examples
% > lit=haslight(gca)
% lit =
%      1

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

lit = false;
ch = get(axish,'Children');

for k = 1:length(ch)
    type = get(ch(k),'Type');
    if strcmp(type,'light')
        lit = true;
        return;
    end
end
%
%end function
