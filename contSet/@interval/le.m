function res = le(obj1,obj2)
% le - Overloads the <= operator;
%    here: Is one interval equal or the subset of another interval?
%
% Syntax:  
%    res = le(obj1,obj2)
%
% Inputs:
%    obj1 - interval object
%    obj2 - another interval object
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

% Author:       Matthias Althoff
% Written:      22-July-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%all left borders of obj1 bigger than obj2?
leftResult = all(infimum(obj1) >= infimum(obj2));

%all right borders of obj1 smaller than obj2?
rightResult = all(supremum(obj1) <= supremum(obj2));

%left and right interval test must be true
res = leftResult & rightResult;

%------------- END OF CODE --------------