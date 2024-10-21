function S_out = copy(S)
% copy - copies the contSet object (used for dynamic dispatch)
%
% Syntax:
%    S_out = copy(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    S_out - copied contSet object
%
% Example: 
%    S = zonotope([1;0],[1 0 -1; 0 1 1]);
%    S_out = copy(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:notSupported',...
    'The chosen subclass of contSet does not support a dynamic copy operation.'));

% ------------------------------ END OF CODE ------------------------------
