function P = minus(P,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference
%
% Syntax:  
%    P = minus(P,varargin)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    P - mptPolytope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope/minkDiff

% Author:       Mark Wetzlinger
% Written:      09-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isnumeric(varargin{1})
    % subtrahend is numeric
    P = minkDiff(P,varargin{:});

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class mptPolytope except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

%------------- END OF CODE --------------
