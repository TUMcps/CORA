function Int = generateRandom(varargin)
% generateRandom - Generates a random interval
%
% Syntax:  
%    Int = generateRandom(varargin)
%
% Inputs:
%    dim - (optional) dimension
%
% Outputs:
%    Int - random interval
%
% Example: 
%    Int = interval.generateRandom();
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 1
    dim_x = varargin{1};
else
    dim_low = 1;
    dim_up  = 10;
    dim_x = dim_low + floor(rand(1) * (dim_up - dim_low + 1));
end

% ranges
range_low = -10;
middle = 0;
range_up = 10;

% instantiate interval
Int = interval(range_low + rand(dim_x,1)*(middle - range_low),...
        middle + rand(dim_x,1)*(range_up - middle));

%------------- END OF CODE --------------