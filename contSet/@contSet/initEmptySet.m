function S = initEmptySet(type)
% initEmptySet - instantiates an empty set of a contSet class
%
% Syntax:
%    S = initEmptySet(type)
%
% Inputs:
%    type - contSet class name
%
% Outputs:
%    S - instantiated set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

throw(CORAerror('CORA:specialError', ...
    sprintf(['The function ''contSet.initEmptySet'' is deprecated (since CORA 2024.1.0) and has been replaced by ''contSet.empty''.\n' ...
             'Note that the function ''initEmptySet'' will be removed in a future release.'])));

% ------------------------------ END OF CODE ------------------------------
