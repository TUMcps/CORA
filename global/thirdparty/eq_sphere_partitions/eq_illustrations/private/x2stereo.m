function result = x2stereo(x)
%X2STEREO Stereographic projection of Euclidean points
%
% result = x2stereo(x);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(x,1)-1;
mask = (x(dim+1,:) == 1);
scale = ones(dim,1)*(1-x(dim+1,~mask));
result = x(1:dim,~mask)./scale;
% end function
