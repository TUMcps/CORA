function res = isscalar(I)
% isscalar - check if an interval is one-dimensional
%
% Syntax:
%    res = isscalar(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I = interval(-1,2);
%    res = isscalar(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       25-June-2015
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check of infimum sufficient
res = isscalar(I.inf);

% ------------------------------ END OF CODE ------------------------------
