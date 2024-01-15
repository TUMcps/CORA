function S = empty(n)
% empty - instantiates an empty set
%
% Syntax:
%    S = empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:notSupported',...
    'The chosen subclass of contSet does not support an empty set instantiation.'));

% ------------------------------ END OF CODE ------------------------------
