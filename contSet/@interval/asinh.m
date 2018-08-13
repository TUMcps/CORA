function res = asinh(intVal)
% asinh - Overloaded 'asinh()' operator for intervals
%
% Syntax:  
%    res = asinh(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      12-February-2016
% Last update:  21-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% x_ is x infimum, x-- is x supremum

% [asinh(x_), asinh(x--)].

%------------- BEGIN CODE --------------


res = interval();

res.inf = asinh(intVal.inf);
res.sup = asinh(intVal.sup);
        


% matrix case
% rand(100, 100)
%time_CORA =
%   0.001107971243036
%time_INTLAB =
%   0.010044490163099


%------------- END OF CODE --------------