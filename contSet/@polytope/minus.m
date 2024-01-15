function P_out = minus(P,varargin)
% minus - dummy function to alert users of the difference in meaning
%    between 'minus' for range bounding and 'minkDiff' for the Minkowski
%    difference; for numerical vectors as subtrahends, these operations are
%    equivalent, so we compute it here nonetheless
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
%    P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
%    v = [1;0];
%    P_ = P - v;
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
    % copy properties
    P_out.bounded.val = P.bounded.val;
    P_out.emptySet.val = P.emptySet.val;
    P_out.fullDim.val = P.fullDim.val;
    P_out.minVRep.val = P.minVRep.val;
    P_out.minHRep.val = P.minHRep.val;

else
    % throw error
    throw(CORAerror('CORA:notSupported',...
        ['The function ''minus'' is not implemented for the class polytope except for vectors as a subtrahend.\n', ...
        'If you require to compute the Minkowski difference, use ''minkDiff'' instead.']));
end

% ------------------------------ END OF CODE ------------------------------
