function res = hashFunction( monomials )
% hashFunction - adds the sum of all columns as the first column
%
% Syntax:  
%    value = hashFunction( monomials )
%
% Inputs:
%    monomials - vector with the multivariate monomials of one terms 
%               (i.e [2 1] for x.^2 * y)
%
% Outputs:
%    value - resulting rank to the term (i.e. 2D-case with max order 11:
%               [2 1] -> 030201
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm
%
% Author:       Niklas Kochdumper, Dmitry Grebenyuk
% Written:      14-June-2017
%               02-December-2017 (DG) New rank evaluation
% Last update:  ---  
% Last revision:---

%------------- BEGIN CODE -------------
 
    res = [sum(monomials, 2), monomials];

end

%------------ END OF CODE ------------