function aff = affine(I,varargin)
% affine - conversion to affine objects
%
% Syntax:
%    I = affine(I)
%
% Inputs:
%    I - interval object
%    name - a cell containing a name of a variable
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models 
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of Taylor models
%    eps - precision for the selected optimization method (opt_method = 'bnb', 
%          opt_method = 'bnbAdv', and opt_method = 'linQuad')
%    tolerance - monomials with coefficients smaller than this value are
%                moved to the remainder
%
% Outputs:
%    aff - affine object
%
% Example:
%    I = interval(0,pi/2);
%    aff = affine(I,'a','int',[],1e-8);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: affine/affine

% Authors:       Mark Wetzlinger
% Written:       23-December-2023 (MW, moved from affine constructor)
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

aff = affine(I.infimum,I.supremum,varargin{:});

% ------------------------------ END OF CODE ------------------------------
