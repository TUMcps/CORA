function res = eq(I1,I2)
% eq - Overloads the == operator; here: Are both intervals equal?
%
% Syntax:  
%    res = eq(I1,I2)
%
% Inputs:
%    I1 - interval object
%    I2 - interval object
%
% Outputs:
%    res - true/false whether intervals are equal
%
% Example: 
%    I1 = interval([1; -1], [2; 1]);
%    I2 = interval([1; -2], [2; 2]);
%    I1 == I2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      05-August-2016 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check equality of infima
leftResult = all(all(infimum(I1) == infimum(I2)));

% check equalify of suprema
rightResult = all(all(supremum(I1) == supremum(I2)));

% combine checks
res = leftResult & rightResult;

%------------- END OF CODE --------------