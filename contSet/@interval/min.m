function res = min(I, Y, varargin)
% min - computes the minimum of I and Y
%
% Syntax:
%    res = min(I, Y)
%
% Inputs:
%    I - interval object
%    Y - interval or numeric 
%    varargin - additional parameters for build-in min function
%
% Outputs:
%    res - false
%
% Example: 
%    I = interval([-2;-1],[2;1]);
%    res = min(I, 0)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: min, interval/infimum, interval/max

% Authors:       Tobias Ladner
% Written:       16-December-2022
% Last update:   11-April-2024 (TL, single input)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    % return infimum
    res = I.inf;
    return
end

% inputArgsCheck in interval/max
res = -max(-I, -1 * Y, varargin{:});

% ------------------------------ END OF CODE ------------------------------
