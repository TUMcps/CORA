function n_collars = num_collars(N,c_polar,a_ideal)
%NUM_COLLARS The number of collars between the polar caps
%
% Given N, an ideal angle, and c_polar,
% determine n_collars, the number of collars between the polar caps.
%
%  n_collars = num_collars(N,c_polar,a_ideal);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_collars = zeros(size(N));
enough = (N > 2) & (a_ideal > 0);
n_collars(enough) = max(1,round((pi-2*c_polar(enough))./a_ideal(enough)));
%
% end function
