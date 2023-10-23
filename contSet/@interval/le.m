function res = le(I1,I2)
% le - Overloads the <= operator, checks whether one interval is a subset
%    or equal to another interval
%
% Syntax:
%    res = le(I1,I2)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%
% Outputs:
%    res - Boolean variable: 1 if obj1 is subset or equal to obj2
%
% Example: 
%    I1 = interval([1; -1], [2; 1]);
%    I2 = interval([1; -2], [2; 2]);
%    I1 <= I2
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

%all left borders of obj1 bigger than obj2?
leftResult = all(infimum(I1) >= infimum(I2));

%all right borders of obj1 smaller than obj2?
rightResult = all(supremum(I1) <= supremum(I2));

%left and right interval test must be true
res = leftResult & rightResult;

% ------------------------------ END OF CODE ------------------------------
