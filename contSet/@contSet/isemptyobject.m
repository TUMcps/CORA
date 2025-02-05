function res = isemptyobject(S)
% isemptyobject - checks whether a contSet object contains any information
%    at all; consequently, the set is equivalent to the empty set 
%
% Syntax:
%    res = isemptyobject(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    isemptyobject(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-June-2022
% Last update:   18-August-2022 (MW, extend to class arrays)
%                24-July-2023 (MW, move checks to subclasses, throw error)
% Last revision: 07-February-2025 (TL, removed old, unused aux checks)

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:noops',S));

end

% ------------------------------ END OF CODE ------------------------------
