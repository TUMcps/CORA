function h = generateRandom(varargin)
% generateRandom - Generates a random halfspace
%
% Syntax:  
%    h = generateRandom(varargin)
%
% Inputs:
%    dim - (optional) dimension
%
% Outputs:
%    h - random halfspace
%
% Example: 
%    h1 = halfspace.generateRandom();
%    h2 = halfspace.generateRandom(3);
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
range_up = 10;

% instantiate interval
h = halfspace(range_low + rand(dim_x,1)*(range_up - range_low),...
        range_low + rand(1)*(range_up - range_low));

%------------- END OF CODE --------------