function C = generateRandom(varargin)
% generateRandom - Generates a random capsule
%
% Syntax:  
%    C = generateRandom(varargin)
%
% Inputs:
%    dim - (optional) dimension
%    cen - (optional) center
%
% Outputs:
%    C - random capsule
%
% Example: 
%    C = capsule.generateRandom();
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

dim_low = 1;
dim_up  = 10;
range_low = -10;
range_up  = 10;
totalrange = range_up - range_low;
maxradius = 5;

if nargin == 2
    dim_x = varargin{1};
    cen = varargin{2};
elseif nargin == 1
    dim_x = varargin{1};
    cen = range_low + rand(dim_x,1)*totalrange;
else
    dim_x = dim_low + floor(rand(1) * (dim_up - dim_low + 1));
    cen = range_low + rand(dim_x,1)*totalrange;
end


C = capsule(cen,range_low + rand(dim_x,1)*totalrange,rand(1)*maxradius);


%------------- END OF CODE --------------