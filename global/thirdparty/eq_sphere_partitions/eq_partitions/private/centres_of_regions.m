function points = centres_of_regions(regions)
%CENTRES_OF_REGIONS Centre points of given regions
%
% points = centres_of_regions(regions);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

tol = eps*2^5;
dim = size(regions,1);
N =   size(regions,3);
points = zeros(dim,N);
top = regions(:,1,:);
bot = regions(:,2,:);
zero_bot = abs(bot(1,:)) < tol;
bot(1,zero_bot) = 2*pi;
equal_bot = abs(bot(1,:) - top(1,:)) < tol;
bot(1,equal_bot) = top(1,equal_bot) + 2*pi;
twopi_bot = abs(bot(1,:) - top(1,:) - 2*pi) < tol;
points(1,twopi_bot) = 0;
points(1,~twopi_bot) = mod((bot(1,~twopi_bot) + top(1,~twopi_bot))/2, 2*pi);
for k = 2:dim
   pi_bot =   abs(bot(k,:) - pi) < tol;
   points(k,pi_bot) = pi;
   zero_top = abs(top(k,:)) < tol;
   points(k,zero_top) = 0;
   all_else = ~(zero_top | pi_bot);
   points(k,all_else) = mod((top(k,all_else) + bot(k,all_else))/2, pi);
end
% end function
