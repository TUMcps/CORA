function res = lt(I1,I2)
% lt - Overloads the '<'-operator, checks whether one interval is the
%    subset of another interval
%
% Syntax:
%    res = lt(I1,I2)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([1; -1], [2; 1]);
%    I2 = interval([1; -2], [2; 2]);
%    I1 < I2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       22-July-2016 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% all infima of I1 bigger than I2?
leftResult = all(infimum(I1) > infimum(I2));

% all suprema of I1 smaller than I2?
rightResult = all(supremum(I1) < supremum(I2));

% both tests must be true
res = leftResult & rightResult;

% ------------------------------ END OF CODE ------------------------------
