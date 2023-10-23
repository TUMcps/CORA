function res = infimum(I)
% infimum - returns the infimum of an interval
%
% Syntax:
%    res = infimum(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - infimum of interval
%
% Example: 
%    I = interval([-1 1], [1 2]);
%    res = infimum(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = I.inf;

% ------------------------------ END OF CODE ------------------------------
