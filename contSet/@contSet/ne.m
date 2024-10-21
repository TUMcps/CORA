function res = ne(S1, S2, varargin)
% ne - overloads the '~=' operatort
%
% Syntax:
%    res = ne(S1,S2)
%
% Inputs:
%    S1 - contSet
%    S2 - contSet
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(S1, S2, varargin{:});

end

% ------------------------------ END OF CODE ------------------------------
