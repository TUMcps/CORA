function res = hashFunction(monomials)
% hashFunction - adds the sum of all columns as the first column
%
% Syntax:
%    value = hashFunction(monomials)
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

% Authors:       Niklas Kochdumper, Dmitry Grebenyuk
% Written:       14-June-2017
% Last update:   02-December-2017 (DG, new rank evaluation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
res = [sum(monomials, 2), monomials];

% ------------------------------ END OF CODE ------------------------------
