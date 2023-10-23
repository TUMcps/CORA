function I = uminus(I)
% uminus - Overloaded '-' operator for single operand
%
% Syntax:
%    I = uminus(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    I - interval object
%
% Example: 
%    I = interval([-2;-3],[5;6]);
%    -I
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   21-May-2022 (MW, simpler computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% store value
temp = I.inf;

% rewrite infimum and supremum
I.inf = -I.sup;
I.sup = -temp;

% ------------------------------ END OF CODE ------------------------------
