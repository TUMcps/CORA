function l = length(I)
% length - Overloads the operator that returns the length
%    of the longest array dimension
%
% Syntax:
%    l = length(I)
%
% Inputs:
%    I - interval object 
%
% Outputs:
%    l - length of the largest array dimension.
%
% Example: 
%    I = interval([-1 1], [1 2]);
%    length(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-November-2015 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% length of infimum and supremum equal
l = length(I.inf);

% ------------------------------ END OF CODE ------------------------------
