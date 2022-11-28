function cPZ = minus(cPZ,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference
%
% Syntax:  
%    cPZ = minus(cPZ,varargin)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    cPZ - conPolyZono object
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
% See also: conPolyZono/minkDiff

% Author:       Mark Wetzlinger
% Written:      09-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isnumeric(varargin{1})
    % subtrahend is numeric
    cPZ = minkDiff(cPZ,varargin{:});

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class conPolyZono except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

%------------- END OF CODE --------------
