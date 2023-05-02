function S = minkDiff(S1,S2,varargin)
% minkDiff - Minkowski difference
%
% Syntax:  
%    S = minkDiff(S1,S2)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object, numerical vector
%    method - (optional) method
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minkDiff

% Author:       Mark Wetzlinger
% Written:      02-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% throw error if minkDiff is not implemented by subclass
throw(CORAerror('CORA:noops',S1,S2));

%------------ END OF CODE ------------