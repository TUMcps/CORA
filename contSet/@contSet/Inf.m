function S = Inf(n)
% Inf - instantiates a fullspace set
%
% Syntax:
%    S = Inf(n)
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
    'The chosen subclass of contSet does not support a fullspace set instantiation.'));

% ------------------------------ END OF CODE ------------------------------
