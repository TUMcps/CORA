function res = jsonencode(I, varargin)
% jsonencode - Converts an interval object to a JSON string
%
% Syntax:
%    res = jsonencode(I)
%
% Inputs:
%    I - interval object
%    varargin - additional parameters for built-in jsonencode function
%
% Outputs:
%    res - JSON string representation of the interval
%
% Example: 
%    I = interval([1;2],[3;4]);
%    res = jsonencode(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: jsonencode

% Authors:       Sebastian Mair
% Written:       20-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res.inf = num2str(I.inf, '%.4f');
res.sup = num2str(I.sup, '%.4f');
res = jsonencode(res);

% ------------------------------ END OF CODE ------------------------------
    