function S = origin(n)
% origin - instantiates a set representing only the origin
%
% Syntax:
%    S = origin(n)
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
% Written:       21-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({{n,'att','numeric',{'scalar','positive','integer'}}});
% is overridden in subclass if implemented; throw error
throw(CORAerror('CORA:notSupported',...
    'The chosen subclass of contSet does not support representing only the origin.'));

% ------------------------------ END OF CODE ------------------------------
