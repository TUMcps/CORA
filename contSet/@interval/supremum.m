function res = supremum(I)
% supremum - returns the supremum of an interval
%
% Syntax:
%    res = supremum(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - numerical value
%
% Example: 
%    I = interval([-1;1],[1;2]);
%    res = supremum(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = I.sup;

% ------------------------------ END OF CODE ------------------------------
