function P_out = minus(P,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference
%
% Syntax:
%    P_out = minus(P,varargin)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    P_out - polytope object
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
% See also: zonotope/minkDiff

% Authors:       Mark Wetzlinger
% Written:       09-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isnumeric(varargin{1})
    % subtrahend is numeric
    P_out = minkDiff(P,varargin{:});

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class polytope except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

% ------------------------------ END OF CODE ------------------------------
