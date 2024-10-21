function S = enlarge(S, factor)
% enlarge - enlarges the set by the given factor without changing its center
%
% Syntax:
%    S = enlarge(S,factor)
%
% Inputs:
%    S - contSet object
%    factor - column vector of factors for the enlargement of each dimension
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,2)
inputArgsCheck({ ...
    {S,'att','contSet'}, ...
    {factor,'att','numeric','column'}, ...
})

% shift to origin
c = center(S);
S = S - c;

% enlarge set
S = S .* factor;

% shift back to center
S = S + c;

% ------------------------------ END OF CODE ------------------------------
