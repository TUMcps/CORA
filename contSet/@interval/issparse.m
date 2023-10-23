function res = issparse(I)
% issparse - checks if an interval is sparse
%
% Syntax:
%    res = issparse(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I = interval(speye(2),2*speye(2));
%    res = issparse(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%return result
res = issparse(I.inf) || issparse(I.sup);

% ------------------------------ END OF CODE ------------------------------
