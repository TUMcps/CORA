function cZ = minus(cZ,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference
%
% Syntax:  
%    cZ = minus(cZ,varargin)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    ---
%
% References:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/minkDiff

% Author:       Mark Wetzlinger
% Written:      09-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isnumeric(varargin{1})
    % subtrahend is numeric
    cZ = minkDiff(cZ,varargin{:});

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class conZonotope except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

%------------- END OF CODE --------------
